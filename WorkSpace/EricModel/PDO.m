addpath(genpath('../../src'));
addpath('tmpPropensityFunctions');

loadPrevious = false;
savedWorkspace = 'workspaceJuly24';

load('workspaceJuly24.mat');

%%
%% Additional Codes for FIM and PDO analyses (ALEX Paper)
%%  STEP XX -- Explore Optimal Designs Using FIM Add FIM to MH UQ Plots 
%%    STEP XX.A. -- FIM Analyses
%%      STEP XX.A.1. -- Plot UQ from FIM compared to MH

figNew = figure;

GRDusp1_log10PriorStd = log10PriorStd(1:13);

fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
    diag(GRDusp1_log10PriorStd.^2));
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

%%      STEP XX.A.2. -- Find optimal experiment design (same number of cells)
nTotal = sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOpt = ModelGRDusp100nM.optimizeCellCounts(fimResults,nTotal,'tr[1:4]'); % min. inv determinant <x^{-1}> (all other parameters are known and fixed)
nCellsOptAvail = min(nCellsOpt,ModelGRDusp100nM.dataSet.nCells')
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));
fimOptAvail = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptAvail,diag(GRDusp1_log10PriorStd.^2));
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


%%    STEP XX.B. -- PDO Calculations
%%      STEP XX.B.1. -- Calibrate PDO from Multi-Modal Experimental Data
% Calibration the PDO from empirical data. Here, the number of spots has
% been measured using different assays in data columns 'nTotal' for the
% 'true' data set and in the columns 'nSpots0' for a different label or
% 'intens1' for the integrated intensity.  We calibrate two different PDOs
% for this case. In both cases, we assume an 'AffinePoiss' PDO where the
% obervation probability is a Poisson distribution where the mean value is
% affine linearly related to the true value: P(y|x) = Poiss(a0 + a1*x);
ModelPDOSpots = ModelGRDusp100nM.calibratePDO('../../ExampleData/pdoCalibrationData.csv',...
    {'rna'},{'nTotal'},{'nSpots0'},'AffinePoiss',true);

%%      STEP XX.B.2. -- Calibrate PDO from Eric's DUSP1 Intensity Data.
ModelPDOIntensEric = ModelGRDusp100nM;
ModelPDOIntensEric = ModelPDOIntensEric.calibratePDO('EricData/pdoCalibrationData_EricIntensity_ConcHigh.csv',...
    {'rna'},{'RNA_DUSP1_nuc'},{'Nuc_DUSP1_avg_int_tot'},'AffinePoiss',true,[1,230,0.5]);

%%    STEP XX.C. -- FIM + PDO Analyses
%%      STEP XX.C.1. -- Analyze FIM with PDO for MCP/smFISH
fimsPDOSpot = ModelPDOSpots.computeFIM(sensSoln.sens,'log');
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));

nCellsOptPDOspots = ModelPDOSpots.optimizeCellCounts(fimsPDOSpot,nTotal,'tr[1:4]');


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
%%      STEP XX.C.2. -- Analyze FIM with PDO for Intensity only
fimsPDOIntens = ModelPDOIntensEric.computeFIM(sensSoln.sens,'log');
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));

nCellsOptPDOintens = ModelPDOSpots.optimizeCellCounts(fimsPDOIntens,nTotal,'tr[1:4]');

fimPDOIntensAvail = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptAvail,diag(GRDusp1_log10PriorStd.^2));
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

%%      STEP XX.C.3. -- Analyze FIM with PDO for Intensity only more cells
nTimes = 3.71;
fimPDOIntens2x = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt*nTimes,diag(GRDusp1_log10PriorStd.^2));
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


%%    STEP XX.D. -- Validation of FIM Predictions
%%      STEP XX.D.1. -- Fit to FIM selected time points
ModelGRDusp100nM_FIMDesign = ModelGRDusp100nM;
% Set the fitting routine to only consider the time points selected by the
% FIM analysis:
ModelGRDusp100nM_FIMDesign.fittingOptions.timesToFit = nCellsOpt>0;
% Refit the model, but now with only those time points.
DUSP1parsFIMDesign = DUSP1pars;
for i = 1:1
    DUSP1parsFIMDesign = ModelGRDusp100nM_FIMDesign.maximizeLikelihood(...
        DUSP1parsFIMDesign, fitOptions);
    ModelGRDusp100nM_FIMDesign.parameters(1:4,2) = num2cell(DUSP1parsFIMDesign);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign') 
end

%%      STEP XX.D.2. -- Plot fit and predicitons using FIM suggested conditions.
showCases = [1,1,1,1];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsFIMDesign,Dusp1FitCases,showCases)

%%      STEP XX.D.3. -- Fit the DISTORTED intensity measurements at ALL times.
% Here, I will load the intensity data for the nuclear DUSP1.  then I will
% attempt to identify the model from just that data and at just the times
% selected by the FIM.
ModelPDOIntensEric.dataSet = [];
ModelPDOIntensEric = ModelPDOIntensEric.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','Nuc_DUSP1_avg_int_tot'},...
        {'Dex_Conc','100'}); 

%load('EricModelDusp1_MMDex','DUSP1parsIntensity') 

DUSP1parsIntensity = DUSP1parsFIMDesign;

ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1parsIntensity);

for i = 1:fitIters
    DUSP1parsIntensity = ModelPDOIntensEric.maximizeLikelihood(...
        DUSP1parsIntensity, fitOptions);
    ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1parsIntensity);
    %save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsIntensity','DUSP1parsFIMDesignIntensity') 
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsIntensity') 
end
%%        STEP XX.D.3.a. -- Plot predictions when fit to distorted data at ALL times.
showCases = [1,1,1,1];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsIntensity,Dusp1FitCases,showCases)

%%      STEP XX.D.4. -- Fit the DISTORTED intensity measurements at FIM selected times.
% Here, I will load the intensity data for the nuclear DUSP1.  then I will
% attempt to identify the model from just that data and at just the times
% selected by the FIM.
ModelPDOIntensEricFIM = ModelPDOIntensEric;
ModelPDOIntensEricFIM.dataSet = [];
ModelPDOIntensEricFIM = ModelPDOIntensEricFIM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','Nuc_DUSP1_avg_int_tot'},...
        {'Dex_Conc','100'}); 

% Set fitting routine only to consider the time points selected by the FIM.
ModelPDOIntensEricFIM.fittingOptions.timesToFit = nCellsOpt>0;

% Refit the model, but now with only those time points.
%load('EricModelDusp1_MMDex','DUSP1parsFIMDesignIntensity') 

DUSP1parsFIMDesignIntensity = DUSP1parsIntensity;

ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesignIntensity);

for i = 1:5
    DUSP1parsFIMDesignIntensity = ModelPDOIntensEricFIM.maximizeLikelihood(...
        DUSP1parsFIMDesignIntensity, fitOptions);
    ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesignIntensity);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsFIMDesignIntensity') 
end
%%        STEP XX.D.5.a. -- Plot predictions when fit to distorted data at FIM times.
showCases = [1,1,1,1];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsFIMDesignIntensity,Dusp1FitCases,showCases)

%%    STEP XX.E. -- Plot Information vs. Expt Design, PDO, and Number of Cells
fimlog10PriorStd = 2*ones(1,13);
fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,diag(fimlog10PriorStd.^2));
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(fimlog10PriorStd.^2));
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(fimlog10PriorStd.^2));
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(fimlog10PriorStd.^2));

barWithOriginalNumber = [det(fimOrig{1}(1:4,1:4)),det(fimOpt{1}(1:4,1:4)),det(fimPDOSpots{1}(1:4,1:4)),det(fimPDOIntens{1}(1:4,1:4))];

nCellVec = logspace(3,5,20);
nCellsOrigRat = ModelGRDusp100nM.dataSet.nCells/sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOptRat = nCellsOpt/sum(nCellsOpt);

cols = {'k','c','b','r'};
fimDetVsNumberMAt=[];
for i = 1:length(nCellVec)
    fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOrigRat*nCellVec(i),diag(fimlog10PriorStd.^2));
    fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptRat*nCellVec(i),diag(fimlog10PriorStd.^2));
    fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOptRat*nCellVec(i),diag(fimlog10PriorStd.^2));
    fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptRat*nCellVec(i),diag(fimlog10PriorStd.^2));
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


%% Extra Functions
function makeGRPlots(combinedModel,GRpars)
combinedGRModel = combinedModel.updateModels(GRpars,false);
nMods = length(combinedGRModel.SSITModels);
ModelGroup = cell(nMods,1);
for i=1:nMods
    %  Update parameters in original models.
    ModelGroup{i} = combinedGRModel.SSITModels{i};
    ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
    ModelGroup{i}.makeFitPlot([],1,[],true,'STD')
end
end

function plotODEresults(extendedMod,soln,modeWithGRData,fignum)
arguments
    extendedMod
    soln
    modeWithGRData
    fignum = 1;
end
figure(fignum); clf;
% Plot GR levels vs. Time
subplot(2,1,1)
plot(extendedMod.tSpan,soln.ode(:,3:4),'--','LineWidth',2);hold on
plot(modeWithGRData.dataSet.times,modeWithGRData.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
legend(extendedMod.species(3:4))
set(gca,'xlim',[-10,200],'ylim',[0,12],'fontsize',16)
ylabel('GR Concentrations (UA)')
legend({'Cyt-Model','Nuc-Model','Cyt-Data','Nuc-Data'})
title('GR')

% Plot DUSP1 levels vs. Time
subplot(2,1,2)
plot(extendedMod.tSpan,soln.ode(:,5:6),'--','LineWidth',2);hold on
plot(extendedMod.dataSet.times,extendedMod.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
legend(extendedMod.species(3:4))
set(gca,'xlim',[-10,200],'ylim',[0,160],'fontsize',16)
ylabel('GR Concentrations (UA)')
legend({'Nuc-Model','Cyt-Model','Nuc-Data','Cyt-Data'})
title('DUSP1')
end

function makeCytDistPlots(ssaSoln_100,extendedMod,fignum,timeIndsMod,timeIndsDat,speciesIndMod,speciesIndDat)
arguments
    ssaSoln_100
    extendedMod
    fignum = 1;
    timeIndsMod = [];
    timeIndsDat = [];
    speciesIndMod = 1;
    speciesIndDat = 1;
end
figure(fignum); clf;
if isempty(timeIndsMod)
    timeIndsMod = [1:length(ssaSoln_100.T_array)];
end
if isempty(timeIndsDat)
    timeIndsDat = [1:length(extendedMod.dataSet.times)];
end

if length(timeIndsMod)~=length(timeIndsDat)
    error('Length of data and model time points must be equal')
end

Nrows = ceil(sqrt(length(timeIndsMod)));
for i = 1:length(timeIndsMod)
    subplot(Nrows,Nrows,i)

    % Add SSA to histogram plot
    M = squeeze(ssaSoln_100.trajs(speciesIndMod,timeIndsMod(i),:));
    H = histogram(M,'Normalization','pdf');
    hold on

    % Add data to histogram plot
    dMat = double(extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(timeIndsDat(i),:,:));
    N = sum(dMat,"all");
    if speciesIndDat==1
        PD = [0;sum(dMat,2)/N];
    else
        PD = [0;sum(dMat,1)'/N];
    end
    
    binEdges = round(H.BinEdges);
    nBins = length(binEdges)-1;
    PDbinned = zeros(nBins,1);
    binwidth = binEdges(2)-binEdges(1);
    for j = 1:nBins
        PDbinned(j) = sum(PD(binEdges(j)+1:binEdges(j+1)));
    end

    PDbinned = [PDbinned;1-sum(PDbinned)];
    stairs(binEdges,PDbinned/binwidth,'linewidth',2)

    set(gca,'FontSize',15,'ylim',[0,0.03])
    
end
end

function [TSLikelihood,modelResults,dataResults] = computeTSlikelihood(x,ModelTS,allData,dex,outlierThresh,reportErrors)
arguments
    x
    ModelTS
    allData
    dex = 100
    outlierThresh = 50;  
    reportErrors = false
    % The outlierThresh input sets an outlier threshold for the nascent RNA
    % at the TS.  Data with observed counts more than this value will be
    % trincated to this value.  Model predictions for counts higher than
    % this value will also be integrated into a single bin "more than
    % outlierThresh".
end

tau_elong = x(end);  % Time for elongation and processing of mRNA (assume deterministic)

% The rest of the parameters go into the model
x = x(1:end-1);
ModelTS.parameters(ModelTS.fittingOptions.modelVarsToFit,2) = num2cell(x);
ModelTS.parameters(13,:) = {'Dex0',dex};

% set upper bound on FSP solution
ModelTS.fspOptions.bounds(6) = outlierThresh*1.5;

% turn off FSP expansion routine
ModelTS.fspOptions.fspTol = inf;

% Create past model only used to get distribution of states for when
% elongation would begin.
ModelTSPast = ModelTS;
ModelTSPast.parameters(3,:) = {'kr',0};  % Turn off transcription.
ModelTSPast.initialCondition(5) = 0;

% Shift solution time points back by tau_elong.
ModelTSPast.tSpan = ModelTS.tSpan - tau_elong;
ModelTSPast.tSpan = [ModelTSPast.tSpan(1)-10,ModelTSPast.tSpan];
ModelTSPast.initialTime = ModelTSPast.tSpan(1);

% Solve for distributions at t - tau_elong.
TSPastSoln = ModelTSPast.solve;

% For the present model we start one elongation period in the past and then 
% integrate forward in time.
ModelTSPresent = ModelTS;
ModelTSPresent.parameters(4,:) = {'knuc2cyt',0};
ModelTSPresent.parameters(14,:) = {'knucdeg',0};

% We will NOT use steady state initial conditions -- rather, we use the
% computed distributions from previous time step.
ModelTSPresent.fspOptions.initApproxSS = false;
% We need to solve once for each time point, where the initial condition is
% defined by the past model.
NT = length(TSPastSoln.fsp)-1;
fsp = cell(1,NT);
for iT = 1:NT
    tmpModelTSPresent = ModelTSPresent.setICfromFspVector(TSPastSoln.stateSpace,TSPastSoln.fsp{iT+1});
    tmpModelTSPresent.tSpan = ModelTSPast.tSpan(iT+1)+[0,tau_elong/2,tau_elong];
    tmpModelTSPresent.initialTime = tmpModelTSPresent.tSpan(1);
    tsPresentSoln = tmpModelTSPresent.solve;
    fsp(iT) = tsPresentSoln.fsp(end); 
end

% Summarize the model results for plotting and for likelihood function
% calculations.
modelResults.meanNascentMod = zeros(NT,1);
modelResults.fracTSMod = zeros(NT,1);
TSthresh = 4;
for iT = NT:-1:1
    P=double(fsp{iT}.p.sumOver([1,2]).data);
    P = max(0,P);
    if sum(P)>1
        P = P/sum(P);
    end
    P(1:TSthresh) = [sum(P(1:TSthresh));zeros(TSthresh-1,1)]; 
    modelResults.Psaved{iT} = P;
end

% Summarize the data for plotting and for likelihood function
% calculations.
for iT = NT:-1:1
    time = ModelTS.tSpan(iT);
    redData = allData(allData.Dex_Conc==dex&allData.Time_index==time,:);
    TS_counts0 = redData.DUSP1_ts_size_0;
    TS_counts1 = redData.DUSP1_ts_size_1;    
    TS_counts0(isnan(TS_counts0)) = 0;
    TS_counts1(isnan(TS_counts1)) = 0;
    dataResults.TS_counts{iT} = TS_counts0+TS_counts1; 
    dataResults.meanNascentDat(iT) = mean(dataResults.TS_counts{iT});
    dataResults.fracTSDat(iT) = sum(dataResults.TS_counts{iT}>=TSthresh)/length(dataResults.TS_counts{iT});

    if reportErrors
        reps = unique(redData.Replica);
        frac = [];
        mn = [];
        for j = 1:length(reps)
            AllDatRedReps = redData(strcmp(redData.Replica,reps{j}),:);
            frac(j) = sum(AllDatRedReps.DUSP1_ts_size_0>0)/size(AllDatRedReps,1);
            mn(j) = (sum(AllDatRedReps.DUSP1_ts_size_0(AllDatRedReps.DUSP1_ts_size_0>0))+...
                sum(AllDatRedReps.DUSP1_ts_size_1(AllDatRedReps.DUSP1_ts_size_1>0)))/size(AllDatRedReps,1);
        end
        dataResults.fracTSDatstd(iT) = std(frac);
        dataResults.meanNascentDatstd(iT) = std(mn);
    end
end


% Adjust to account for partial transcripts.
nRNAmax = floor(ModelTS.fspOptions.bounds(6)+1);
PDO_partial_transcripts = zeros(nRNAmax,nRNAmax);
for i = 1:nRNAmax
    PDO_partial_transcripts(1:i,i) = binopdf(0:i-1,i-1,0.5);
end

% Calculate the likelihood function from the TS data.
TSLikelihood=0;
for iT=1:NT
    modelResults.Psaved{iT} = max(1e-10,modelResults.Psaved{iT});
    
    % Remove outliers from data and set to threshold.
    countsRNA = dataResults.TS_counts{iT};
    countsRNA(countsRNA>outlierThresh) = outlierThresh;
   
    % Apply PDO to account for partial transcripts.
    try
        modelResults.Psaved{iT} = PDO_partial_transcripts*modelResults.Psaved{iT};
    catch
        1+1
    end
    % Bin outlier predictions for model
    P = modelResults.Psaved{iT};

    P(outlierThresh+1) = sum(P(outlierThresh+1:end));
    P = P(1:outlierThresh+1);

    modelResults.meanNascentMod(iT) = [TSthresh:length(P)-1]*P(TSthresh+1:end);
    modelResults.fracTSMod(iT) = sum(P(TSthresh+1:end));

    TSLikelihood = TSLikelihood + sum(log(P(countsRNA+1)));
end
end

function makeTSPlots(modelResults,dataResults,ModelTS,TSthresh,figNum)
arguments
    modelResults
    dataResults
    ModelTS
    TSthresh = 4;
    figNum = [700,701];
end
figure(figNum(1)); clf
NT = length(modelResults.meanNascentMod);
for iT = 1:NT
    P=modelResults.Psaved{iT};
   
    % Lump predictions less than detection limit at zero.
    P(1:TSthresh) = sum(P(1:TSthresh))/4;  
    % P(1:TSthresh) = [sum(P(1:TSthresh));zeros(TSthresh-1,1)];   
    
    subplot(ceil(NT/3),3,iT)
    % stairs([0:length(P)-1],cumsum(P),'linewidth',3)
    stairs([0:length(P)-1],(P),'linewidth',3)
    set(gca,'fontsize',16,'ylim',[0,1],'xlim',[0,50])
    set(gca,'fontsize',16,'ylim',[0,0.02],'xlim',[0,50])
    if iT>9
        xlabel('Number nascent RNA')
    end
    if mod(iT-1,3)==0
        ylabel('Probability')
    end   
end

figure(figNum(2)); clf
subplot(2,1,1)
plot(ModelTS.tSpan,modelResults.meanNascentMod,'linewidth',3);
set(gca,'fontsize',16,'ylim',[0,3])
ylabel('Mean Nascent RNA')

subplot(2,1,2)
plot(ModelTS.tSpan,modelResults.fracTSMod,'linewidth',3);
set(gca,'fontsize',16,'ylim',[0,0.5])
ylabel('Fraction with TS')
xlabel('time (min)')
%
%    STEP 4.C.2. -- Add data to plots
figure(figNum(1))
for iT = 1:NT
    subplot(ceil(NT/3),3,iT); hold on
    histogram(dataResults.TS_counts{iT},[0,TSthresh:52],'Normalization','pdf','DisplayStyle','bar','LineWidth',1);
    % histogram(dataResults.TS_counts{iT},[0,TSthresh:52],'Normalization','cdf','DisplayStyle','stairs','LineWidth',3);
    % nCount(iT) = length(dataResults.TS_counts{iT});
    % semCount(iT) = std(dataResults.TS_counts{iT})/sqrt(nCount(iT)-1);
end

figure(figNum(2))
subplot(2,1,1)
hold on
% plot(ModelTS.tSpan,dataResults.meanNascentDat,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
errorbar(ModelTS.tSpan,dataResults.meanNascentDat,dataResults.meanNascentDatstd,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');

subplot(2,1,2)
hold on
% stdFrac = sqrt(nCount.*dataResults.fracTSDat.*(1-dataResults.fracTSDat))./nCount;
% plot(ModelTS.tSpan,dataResults.fracTSDat,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
errorbar(ModelTS.tSpan,dataResults.fracTSDat,dataResults.fracTSDatstd,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
end