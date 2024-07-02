%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations
close all 
clear all
addpath(genpath('../../src'));

%% Create Base Model for GR Only
% Here, I set up the model for the GR translocation dynamics.
ModelGR = SSIT;
ModelGR.species = {'cytGR';'nucGR'};
ModelGR.initialCondition = [20;1];

ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';...
    'kg1';'gg1*cytGR';'gg2*nucGR'};
ModelGR.stoichiometry = [-1,1,1,-1,0;...
                         1,-1,0,0,-1];

ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
    'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',14e-5;...
    'gg1',1e-5;'gg2',1e-6;'MDex',5;'Dex0',100});

log10PriorMean = [-1 -1 0 -2,...
    -1 -3 -2 -1 -2 -2 -2 0.5, 2];
log10PriorStd = 2*ones(1,13);

ModelGR.fspOptions.initApproxSS = true;

ModelGR.fittingOptions.modelVarsToFit = (5:12);
ModelGR.fittingOptions.logPrior = @(x)-sum((x-log10PriorMean(5:12)).^2./(2*log10PriorStd(5:12).^2));

ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
ModelGR = ModelGR.formPropensitiesGeneral('EricRonModGR');
[FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve;
[FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve(FSPGrSoln.stateSpace);

%%    Load previously fit parameter values (optional)
load('EricModel_MMDex','GRpars')
ModelGR.parameters(5:12,2) = num2cell([GRpars]);

%%    Associate GR Data with Different Instances of Model (10,100nm Dex)
GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
    '10','10',102,'GR Fit (10nM Dex)';...
    '100','100',103,'GR Fit (100nM Dex)'};

ModelGRparameterMap = cell(1,size(GRfitCases,1));
ModelGRfit = cell(1,size(GRfitCases,1));
ModelGRODEfit = cell(1,size(GRfitCases,1));
for i=1:3
    ModelGRfit{i} = ModelGR.loadData("EricDataJan23_2024/Gated_dataframe_Ron_030624_NormalizedGR_bins.csv",...
        {'nucGR','normgrnuc';'cytGR','normgrcyt'},...
        {'Condition','GR_timesweep';'Dex_Conc',GRfitCases{i,2}});
    ModelGRfit{i}.parameters{13,2} = str2num(GRfitCases{i,1});    
    ModelGRparameterMap(i) = {(1:8)};
end

%%    Make Guesses for the FSP bounds
% This is sometimes necessary when using an uninduced steady state as the
% initial condition. You need to guess a reasonalbe statespace or the
% computation of the SS can be inaccurate.
for i = 1:3
    boundGuesses{i} = [0;0;30;30];
end

%%    Combine all three GR models and fit using a single parameter set.
for jj = 1:5
    fitOptions = optimset('Display','iter','MaxIter',300);

    combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap);
    combinedGRModel = combinedGRModel.initializeStateSpaces(boundGuesses);
    combinedGRModel = combinedGRModel.updateModels(GRpars,false);
    GRpars = combinedGRModel.maximizeLikelihood(...
        GRpars, fitOptions);
    save('EricModel_MMDex','GRpars') 
end

%% Compute FIM
% combinedGRModel = combinedGRModel.computeFIMs;

%%    Run MH on GR Models.
% MHFitOptions.thin=1;
% MHFitOptions.numberOfSamples=100;
% MHFitOptions.burnIn=0;
% MHFitOptions.progress=true;
% MHFitOptions.numChains = 1;
% MHFitOptions.useFIMforMetHast = true;
% MHFitOptions.saveFile = 'TMPEricMHGR.mat';
% [~,~,MHResultsGR] = combinedGRModel.maximizeLikelihood(...
%     GRpars, MHFitOptions, 'MetropolisHastings');
% delete(MHFitOptions.saveFile)
 
%%    Make Plots of GR Fit Results
combinedGRModel = combinedGRModel.updateModels(GRpars,false);
for i=1:length(ModelGRfit)
    %  Update parameters in original models.
    ModelGRfit{i} = combinedGRModel.SSITModels{i};
    ModelGRfit{i}.tSpan = sort(unique([ModelGRfit{i}.tSpan,linspace(0,180,30)]));
    ModelGRfit{i}.makeFitPlot([],1,[],true,'STD')
end

%% Extend Model to Include DUSP1 Activation, Production, and Degradation
% Copy parameters from the 100nM Dex stim case in GR.
ModelGRDusp = ModelGRfit{3};
ModelGRDusp.species = {'offGene';'onGene';'cytGR';'nucGR';'rna'};
ModelGRDusp.initialCondition = [2;0;24;1;5];
ModelGRDusp.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
    '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR';...
    'kr*onGene';'gr*rna'};
ModelGRDusp.stoichiometry = [-1,1,0,0,0,0,0,0,0;...
                         1,-1,0,0,0,0,0,0,0;...
                         0,0,-1,1,1,-1,0,0,0;...
                         0,0,1,-1,0,0,-1,0,0;...
                         0,0,0,0,0,0,0,1,-1];
ModelGRDusp.useHybrid = true;
ModelGRDusp.hybridOptions.upstreamODEs = {'cytGR','nucGR'};
ModelGRDusp.solutionScheme = 'FSP';
ModelGRDusp.fspOptions.bounds = [0;0;0;2;2;400];
ModelGRDusp.fittingOptions.modelVarsToFit = 1:4;
ModelGRDusp = ModelGRDusp.formPropensitiesGeneral('EricModDusp1');
duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(1:4)).^2./(2*log10PriorStd(1:4).^2));
ModelGRDusp.fittingOptions.logPrior = duspLogPrior;

%%    Load pre-fit parameters into model.
load('EricModelDusp1_MMDex','DUSP1pars')
ModelGRDusp.parameters(1:4,2) = num2cell(DUSP1pars);

%%    Load and Associate with DUSP1 smFISH Data (100nM Dex Only)
Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
for i = 1:size(Dusp1FitCases,1)
    % TODo - Adjust for newly processed data.
    ModelDusp1Fit{i} = ModelGRDusp.loadData('EricDataJan23_2024/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','RNA_DUSP1_nuc'},...
        {'Dex_Conc','100'}); 
    ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};

    ModelDusp1parameterMap{i} = (1:4);
    % Set Dex concentration.
    ModelDusp1Fit{i}.parameters{13,2} = str2num(Dusp1FitCases{i,1});
    ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
end
DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp.fittingOptions.modelVarsToFit,2}];

ModelGRDusp100nM = ModelGRDusp.loadData('EricDataJan23_2024/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','RNA_DUSP1_nuc'},{'Dex_Conc','100'}); 
%%    Fit DUSP1 model(s) with single parameter set.
% Specify prior on parameters.
for i = 1:5
    fitOptions.suppressFSPExpansion = true;
    DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
        DUSP1pars, fitOptions);
    ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
    ModelGRDusp.parameters(1:4,2) = num2cell(DUSP1pars);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars') 
end

%%    Plot predictions for other Dex concentrations.
showCases = [1,0,0,0];
makePlotsDUSP1(ModelDusp1Fit,ModelGRDusp,DUSP1pars,Dusp1FitCases,showCases)
%% Sample uncertainty for Dusp1 Parameters
%%    Compute sensitivity of the fSP solution
ModelGRDusp100nM.solutionScheme = 'fspSens';
sensSoln = ModelGRDusp100nM.solve();
ModelGRDusp100nM.solutionScheme = 'FSP';
%%    Compute FIM
% define which species in model are not observed.
ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};

% compute the FIM
fimResults = ModelGRDusp100nM.computeFIM(sensSoln.sens,'log');
fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
    diag(log10PriorStd.^2));

FIMfree = fimTotal{1}(1:4,1:4);
if min(eig(FIMfree))<1
    disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
    FIMfree = FIMfree + 1*eye(length(FIMfree));
end
covFree = FIMfree^-1;
covFree = 0.5*(covFree+covFree');
%
%%    Run Metropolis Hastings Search
MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
MHFitOptions.thin=1;
MHFitOptions.numberOfSamples=10000;
MHFitOptions.burnIn=0;
MHFitOptions.progress=true;
MHFitOptions.numChains = 1;
MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
[DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
    [], MHFitOptions, 'MetropolisHastings');
delete('TMPEricMHDusp1.mat')
ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);

%%    Plot the MH results
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1]; %
        CH(1).LineWidth = 3;
    end
end
%%    Add FIM to MH UQ Plots
figNew = figure;

fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
    diag(log10PriorStd.^2));
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimTotal],'log',[],figNew)
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];
        CH(2).LineWidth = 3;
    end
end

%% Find optimal experiment design (same number of cells)
nTotal = sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOpt = ModelGRDusp100nM.optimizeCellCounts(fimResults,nTotal,'tr[1:4]');
nCellsOptAvail = min(nCellsOpt,ModelGRDusp100nM.dataSet.nCells')
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(log10PriorStd.^2));
fimOptAvail = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptAvail,diag(log10PriorStd.^2));
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimTotal],'log',[],figNew);
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOptAvail,fimTotal],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];
        CH(2).LineWidth = 3;
        CH(3).Color=[0,1,1];
        CH(3).LineWidth = 3;
    end
end

f = figure;
set(f,'Position',[616   748   412   170])
bar([1:12],ModelGRDusp100nM.dataSet.nCells,0.45)
hold on
bar([1:12]+0.5,nCellsOpt',0.45)
set(gca,'xtick',[1:12]+0.25,'xticklabel',ModelGRDusp100nM.dataSet.times,'fontsize',16,'ylim',[0,7000])
legend('Intuitive Design','Optimal Design')


%%    Calibrate PDO from Multi-Modal Experimental Data
% Calibration the PDO from empirical data. Here, the number of spots has
% been measured using different assays in data columns 'nTotal' for the
% 'true' data set and in the columns 'nSpots0' for a different label or
% 'intens1' for the integrated intensity.  We calibrate two different PDOs
% for this case. In both cases, we assume an 'AffinePoiss' PDO where the
% obervation probability is a Poisson distribution where the mean value is
% affine linearly related to the true value: P(y|x) = Poiss(a0 + a1*x);
ModelPDOSpots = ModelGRDusp100nM.calibratePDO('../../ExampleData/pdoCalibrationData.csv',...
    {'rna'},{'nTotal'},{'nSpots0'},'AffinePoiss',true);

%%    Calibrate PDO from Eric's DUSP1 Intensity Data.
ModelPDOIntensEric = ModelGRDusp100nM;
ModelPDOIntensEric = ModelPDOIntensEric.calibratePDO('EricDataJan23_2024/pdoCalibrationData_EricIntensity_ConcHigh.csv',...
    {'rna'},{'RNA_DUSP1_nuc'},{'Nuc_DUSP1_avg_int_tot'},'AffinePoiss',true,[1,230,0.5]);

%%  With PDO for MCP/smFISH
fimsPDOSpot = ModelPDOSpots.computeFIM(sensSoln.sens,'log');
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(log10PriorStd.^2));
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimPDOSpots,fimTotal],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];
        CH(2).LineWidth = 3;
        CH(4).Color=[0,1,1];
        CH(4).LineWidth = 3;
        CH(3).Color=[0,0,1];
        CH(3).LineWidth = 3;
    end
end
%%  With PDO for Intensity only
fimsPDOIntens = ModelPDOIntensEric.computeFIM(sensSoln.sens,'log');
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(log10PriorStd.^2));
fimPDOIntensAvail = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptAvail,diag(log10PriorStd.^2));
figNew = figure; clf;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimPDOSpots,fimTotal,fimPDOIntens],'log',[],figNew);
% figNew = figure; clf;
% ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimTotal,fimPDOIntensAvail],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(3).Color=[0,0,0];
        CH(3).LineWidth = 3;
        CH(5).Color=[0,1,1];
        CH(5).LineWidth = 3;
        CH(4).Color=[0,0,1];
        CH(4).LineWidth = 3;
        CH(2).Color=[0,1,0];
        CH(2).LineWidth = 3;
    end
end

%%  With PDO for Intensity only 2x number of cells
nTimes = 3.71;
fimPDOIntens2x = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt*nTimes,diag(log10PriorStd.^2));
det(fimOpt{1}(1:4,1:4))/det(fimPDOIntens2x{1}(1:4,1:4))
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimPDOSpots,fimTotal,fimPDOIntens,fimPDOIntens2x],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(4).Color=[0,0,0];
        CH(4).LineWidth = 3;
        CH(6).Color=[0,1,1];
        CH(6).LineWidth = 3;
        CH(5).Color=[0,0,1];
        CH(5).LineWidth = 3;
        CH(3).Color=[0,1,0];
        CH(3).LineWidth = 3;
        CH(2).Color=[1,0,0];
        CH(2).LineWidth = 3;
    end
end



%% Attempt Fit to just the selected time points
ModelGRDusp100nM_FIMDesign = ModelGRDusp100nM;
% Set the fitting routine to only consider the time points selected by the
% FIM analysis:
ModelGRDusp100nM_FIMDesign.fittingOptions.timesToFit = nCellsOpt>0;
% Refit the model, but now with only those time points.
% DUSP1parsFIMDesign = DUSP1pars;
for i = 1:1
    DUSP1parsFIMDesign = ModelGRDusp100nM_FIMDesign.maximizeLikelihood(...
        DUSP1parsFIMDesign, fitOptions);
    ModelGRDusp100nM_FIMDesign.parameters(1:4,2) = num2cell(DUSP1parsFIMDesign);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign') 
end

%%    Let's see how the predictions do with these parameters.
showCases = [1,1,1,0];
makePlotsDUSP1(ModelDusp1Fit,ModelGRDusp,DUSP1parsFIMDesign,Dusp1FitCases,showCases)

%% Attempt to fit the distorted intensity measurements at ALL times.
% Here, I will load the intensity data for the nuclear DUSP1.  then I will
% attempt to identify the model from just that data and at just the times
% selected by the FIM.
ModelPDOIntensEric.dataSet = [];
ModelPDOIntensEric = ModelPDOIntensEric.loadData('EricDataJan23_2024/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','Nuc_DUSP1_avg_int_tot'},...
        {'Dex_Conc','100'}); 
load('EricModelDusp1_MMDex','DUSP1parsIntensity') 

ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1parsIntensity);

for i = 1:1
    DUSP1parsIntensity = ModelPDOIntensEric.maximizeLikelihood(...
        DUSP1parsIntensity, fitOptions);
    ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1parsIntensity);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsIntensity','DUSP1parsFIMDesignIntensity') 
end

%% Attempt to fit the distorted intensity measurements at FIM selected times.
% Here, I will load the intensity data for the nuclear DUSP1.  then I will
% attempt to identify the model from just that data and at just the times
% selected by the FIM.
ModelPDOIntensEricFIM = ModelPDOIntensEric;
ModelPDOIntensEricFIM.dataSet = [];
ModelPDOIntensEricFIM = ModelPDOIntensEricFIM.loadData('EricDataJan23_2024/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','Nuc_DUSP1_avg_int_tot'},...
        {'Dex_Conc','100'}); 
% Set fitting routine only to consider the time points selected by the FIM.
ModelPDOIntensEricFIM.fittingOptions.timesToFit = nCellsOpt>0;
% Refit the model, but now with only those time points.
load('EricModelDusp1_MMDex','DUSP1parsFIMDesignIntensity') 
ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesignIntensity);

for i = 1:5
    DUSP1parsFIMDesignIntensity = ModelPDOIntensEricFIM.maximizeLikelihood(...
        DUSP1parsFIMDesignIntensity, fitOptions);
    ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesignIntensity);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsFIMDesignIntensity') 
end
%%    Let's see how the predictions do with these parameters.
showCases = [1,1,1,0];
makePlotsDUSP1(ModelDusp1Fit,ModelGRDusp,DUSP1parsFIMDesignIntensity,Dusp1FitCases,showCases)

%% Plot Differences in Information.
fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,diag(log10PriorStd.^2));
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(log10PriorStd.^2));
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(log10PriorStd.^2));
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(log10PriorStd.^2));

barWithOriginalNumber = [det(fimTotal{1}(1:4,1:4)),det(fimOpt{1}(1:4,1:4)),det(fimPDOSpots{1}(1:4,1:4)),det(fimPDOIntens{1}(1:4,1:4))];

nCellVec = logspace(3,5,20);
nCellsOrigRat = ModelGRDusp100nM.dataSet.nCells/sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOptRat = nCellsOpt/sum(nCellsOpt);

cols = {'k','c','b','r'};
fimDetVsNumberMAt=[];
for i = 1:length(nCellVec)
    fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOrigRat*nCellVec(i),diag(log10PriorStd.^2));
    fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptRat*nCellVec(i),diag(log10PriorStd.^2));
    fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOptRat*nCellVec(i),diag(log10PriorStd.^2));
    fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptRat*nCellVec(i),diag(log10PriorStd.^2));
    fimDetVsNumberMAt(i,:) = [det(fimOrig{1}(1:4,1:4)),det(fimOpt{1}(1:4,1:4)),det(fimPDOSpots{1}(1:4,1:4)),det(fimPDOIntens{1}(1:4,1:4))];
end

close all
figNew = figure;
% loglog(nCellVec,1./fimDetVsNumberMAt,'linewidth',2);
for i = 1:4
    loglog([10,1e6],1/barWithOriginalNumber(i)*[1,1],[cols{i},'--'],'linewidth',2)
    hold on
    loglog(nCellVec,1./fimDetVsNumberMAt(:,i),cols{i},'linewidth',2)
end
set(gca,'fontsize',15,'ylim',10.^[-12.,-9],'xlim',10.^[3,5])

%% Stochastic Simulation of Full Model
% Form version of full model to run with SSA.  This will allow us to run
% the full model and test to see if it matches to the FSP Model.
SSAModel = ModelGRDusp100nM;
SSAModel.solutionScheme = 'SSA';
SSAModel.useHybrid = false;
SSAModel.ssaOptions.useParalel = true;
SSAModel.inputExpressions = {'IDex','Dex0*exp(-gDex*t)*(t>=0)'};
SSAModel.tSpan = [-100,ModelGRDusp100nM.tSpan];
SSAModel.initialCondition = [2;0;10;3;20];
SSAModel.initialTime = -100;
SSAModel = SSAModel.formPropensitiesGeneral('SSAEric');
%% Plot model
tic
ssaSoln = SSAModel.solve;
toc
% Plot SSA results
SSAModel.makePlot(ssaSoln,'meansAndDevs')

%% Extend model to include cytoplasmic mRNA
SSAModelCytoDusp1 = SSAModel;
SSAModelCytoDusp1.species = {'offGene';'onGene';'cytGR';'nucGR';'rna';'rCyt'};
SSAModelCytoDusp1.initialCondition = [2;0;24;1;5;50];
SSAModelCytoDusp1.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
    '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR';...
    'kr*onGene';'gnuc*rna';...
    'ktr*rna';'gcyt*rCyt'};
SSAModelCytoDusp1.parameters(4,:) = {'gnuc',0.01}; % add parameter for cytoplasmic decay
SSAModelCytoDusp1.parameters(14,:) = {'gcyt',0.005}; % add parameter for cytoplasmic decay
SSAModelCytoDusp1.parameters(15,:) = {'ktr',0.01896}; % add parameter for cytoplasmic decay
SSAModelCytoDusp1.stoichiometry = [-1,1,0,0,0,0,0,0,0,0,0;...
                         1,-1,0,0,0,0,0,0,0,0,0;...
                         0,0,-1,1,1,-1,0,0,0,0,0;...
                         0,0,1,-1,0,0,-1,0,0,0,0;...
                         0,0,0,0,0,0,0,1,-1,-1,0;...
                         0,0,0,0,0,0,0,0,0,1,-1];

SSAModelCytoDusp1 = SSAModelCytoDusp1.formPropensitiesGeneral('SSAEricCyt');
ssaSolnCyt = SSAModelCytoDusp1.solve;

SSAModelCytoDusp1.makePlot(ssaSolnCyt,'meansAndDevs')

%% Make plots of cytoplasmic RNA.
tempCytModel = ModelGRDusp100nM;
tempCytModel = tempCytModel.loadData('EricDataJan23_2024/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','RNA_DUSP1_cyto'},...
        {'Dex_Conc','100'}); 
tempCytModel.makeFitPlot

%% Sandbox for Predicting Other Behaviors
% SBModel = ModelGRDusp;
% SBModel.parameters(13,:) = {'Dex0',100};
% 
% SBModel.tSpan = linspace(0,300,50);
% SBModel.inputExpressions = {'IDex','Dex0*exp(-gDex*t)*(t>0)-Dex0*exp(-gDex*t)*(t>10)'};
% SBModel = SBModel.formPropensitiesGeneral('Pulse');
% [fspSoln] = SBModel.solve;
% SBModel.makePlot(fspSoln,'meansAndDevs',[],[],[4])
% 
% SBModel2 = SBModel;
% SBModel2.tSpan = linspace(0,300,50);
% SBModel2.inputExpressions = {'IDex','Dex0*exp(-gDex*t)*(t>0)'};
% SBModel = SBModel.formPropensitiesGeneral('Step');
% [fspSoln2] = SBModel2.solve;
% SBModel2.makePlot(fspSoln2,'meansAndDevs',[],[],[4])
% 
% %%
% SBModelJoint = ModelGRDusp;
% SBModelJoint.useHybrid = false;
% SBModelJoint.fspOptions.verbose = true;
% SBModelJoint = SBModelJoint.formPropensitiesGeneral('EricGRDusp1Joint');
% SBModelJoint.customConstraintFuns = {'x3+x4','x5/(x4+1)'};
% [fspSoln3] = SBModelJoint.solve;
% SBModelJoint.makePlot(fspSoln3,'joints',[],[])
%% 
