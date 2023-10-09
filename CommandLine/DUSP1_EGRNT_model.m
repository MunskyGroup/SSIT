%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
clear 
close all
clc
addpath('../CommandLine');
addpath('../EricModel/DUSP1_GR_dataframes');
addpath(genpath('../src'));
addpath('../tensor_toolbox-v3.2.1')

%% Model Analysis Options
% Use previous saved results
useSavedModel = false;
useSavedModelMH = false;
useSavedModelSensitivty = false;
useSavedModelMLE = false; 

useManualModelCreation = true; % Create model from scratch
fitModel = true; % Fit model with fminsearch
runMetHast = true; % Calculation of Metropolis Hastings 

calcMLE = false; % MLE calculation for FIM verification
runOptTriptolideExp = false; % optimal triptolide application calc

showFigures = true;

% Save Model, MHresults, sens, and MLE results
% Previous results will not be overwritten
saveMHOutput = true; 
saveModelOutput = true; 
saveSensOutput = true; 
saveMLEOutput = true; 
saveTriptOutput = false;

%% Define SSIT Model (EGRNT)
% Initial Set up for the model
if useManualModelCreation
    EGRNT_Model = SSIT;
    EGRNT_Model.species = {'x1';'x2';'x3'};  % x1: nuclear GR x2: number of on alleles,  x3: mRNA 
    EGRNT_Model.initialCondition = [0;0;0];
    EGRNT_Model.propensityFunctions = {'(kcn0+kcn1*IDex)';'knc*x1';...
                                 'kon*x1*(2-x2)';'koff*x2';'kr*x2';'gr*x3'};
    
    stoich = [1,0,0;   % GR enters nucleus
             -1,0,0;   % GR exist nucleus / degrades
              0,1,0;   % Gene allele switching to on state
              0,-1,0;  % Gene allele switching to off state
              0,0,1;   % mRNA production
              0,0,-1]; % mRNA degredation
    
    EGRNT_Model.stoichiometry = transpose(stoich);
    
    % Input expression for time varying signals
    EGRNT_Model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'}; % Expression for Dex stimulation
    
    % Defining the values of each parameter used
    EGRNT_Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
                         'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
    
    EGRNT_Model.fspOptions.initApproxSS = true;

    EGRNT_Model = EGRNT_Model;
    tpt_array = 0:20:180;
    EGRNT_Model.tSpan = tpt_array; % Define the time points

end



%% Load saved mHast and Sensitivity
if useSavedModel
    j=1;
    while exist(['EGRNT_Model','_v',num2str(j),'.mat'],'file')
        j=j+1;
    end

    if exist(['EGRNT_Model','_v',num2str(j-1),'.mat'],'file')
        EGRNT_Model = load(['EGRNT_Model','_v',num2str(j-1),'.mat']).EGRNT_Model;
    else
        error(['No file named ','EGRNT_Model','_v',num2str(j-1),'.mat found'])
    end

end

if useSavedModelMH
    i=1;
    while exist(['EGRNT_mhast','_v',num2str(i),'.mat'],'file')
        i=i+1;
    end
    if exist(['EGRNT_mhast','_v',num2str(i-1),'.mat'],'file')
        mhResults = load(['EGRNT_mhast','_v',num2str(i-1),'.mat']).mhResults;
    else
        error(['No file named ','EGRNT_mhast','_v',num2str(i),'.mat',' file found'])
    end

end

if useSavedModelSensitivty
    k=1;
    while exist(['EGRNT_sens','_v',num2str(k),'.mat'],'file')
        k=k+1;
    end
    if exist(['EGRNT_sens','_v',num2str(k-1),'.mat'],'file')
        sensSoln = load(['EGRNT_sens','_v',num2str(k-1),'.mat']).sensSoln;
    else
        error(['No file named ','EGRNT_sens','_v',num2str(k-1),'.mat found'])
    end
end

if useSavedModelMLE
    n=1;
    while exist(['EGRNT_MLE','_v',num2str(n),'.mat'],'file')
        n=n+1;
    end
    if exist(['EGRNT_MLE','_v',num2str(n-1),'.mat'],'file')
        EGRNT_MLE = load(['EGRNT_MLE','_v',num2str(n-1),'.mat']).EGRNT_MLE;
        B = load(['B.mat']).B;
    else
        error(['No file named ','EGRNT_MLE','_v',num2str(n),'.mat found'])
    end
end
%% Solve the model using the FSP
EGRNT_Model = EGRNT_Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'}); % Load experimental data set
if fitModel
    EGRNT_Model.fittingOptions.modelVarsToFit = 1:8; % Picking parameters 1-8 for fitting
    
    for i=1:5
        EGRNT_Model.solutionScheme = 'FSP';
        EGRNT_Model.fspOptions.fspTol = 1e-4;
        EGRNT_Model.fspOptions.verbose = 0;
        EGRNT_Model.fspOptions.bounds=[];
        [fspSoln,EGRNT_Model.fspOptions.bounds] = EGRNT_Model.solve;
        EGRNT_Model.fspOptions.bounds
        
        % Load and Fit smFISH Data
        EGRNT_Model.fspOptions.fspTol = inf;
        fitOptions = optimset('Display','iter','MaxIter',500); 
        EGRNT_Model.parameters(EGRNT_Model.fittingOptions.modelVarsToFit,2) =...
          num2cell(EGRNT_Model.maximizeLikelihood(...
          [EGRNT_Model.parameters{EGRNT_Model.fittingOptions.modelVarsToFit,2}],...
          fitOptions));
    end
    
    % figure slide 7 and 8
    if showFigures
        EGRNT_Model.makeFitPlot;
    end
end
%% Metropolis Hastings to Quantify Parameter Uncertainty
if runMetHast
    EGRNT_Model.fittingOptions.modelVarsToFit = 1:4; % Picking parameters 1-4 for fitting
    
    allFitOptions.CovFIMscale = 0.1;% make smaller for higher acceptance rate
    MHOptions = struct('numberOfSamples',15000,'burnin',100,'thin',1,...
      'useFIMforMetHast',true,'suppressFSPExpansion',true);
    [bestParsFound,~,mhResults] = EGRNT_Model.maximizeLikelihood([EGRNT_Model.parameters{EGRNT_Model.fittingOptions.modelVarsToFit,2}]',...
      MHOptions,'MetropolisHastings');
    EGRNT_Model.parameters(EGRNT_Model.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
    
    % figure slide 10
    if showFigures
        EGRNT_Model.plotMHResults(mhResults); 
    end
end

%% Calculate CME Sensitivity to Parameter Variations
EGRNT_Model.sensOptions.solutionMethod = 'finiteDifference';
EGRNT_Model.solutionScheme = 'fspSens'; % Choosing sensitivity solution method
EGRNT_Model.fspOptions.fspTol = 1e-6;
[sensSoln,EGRNT_Model.fspOptions.bounds] = EGRNT_Model.solve;

%% Solve for Fisher Information Matrix at all Time Points
EGRNT_Model.pdoOptions.unobservedSpecies = {'x1','x2'}; % PDO applied for the case where the Gene state is not observed

% FIM calculation taking into account the number of cells measured at each time point
fims = EGRNT_Model.computeFIM(sensSoln.sens);
FIM = EGRNT_Model.evaluateExperiment(fims,EGRNT_Model.dataSet.nCells);

nTotal_dusp = sum(EGRNT_Model.dataSet.nCells);
nCellsOpt_dusp = EGRNT_Model.optimizeCellCounts(fims,nTotal_dusp,'TR[1:4]');
fimOpt_dusp = EGRNT_Model.evaluateExperiment(fims,nCellsOpt_dusp); % Using the FIM to optimize the cell measuremet times

if showFigures
    EGRNT_Model.plotMHResults(mhResults,fimOpt_dusp);
end

%% Maximum likelihood Estimator Verification
if calcMLE
    % Simmulate data with known parameters and re-fit simmulated data to get MLE 
    starttime=-1500;
    Modelssa=EGRNT_Model;
    Modelssa.initialTime = starttime;
    Modelssa.tSpan = [starttime,EGRNT_Model.tSpan];
    Modelssa.ssaOptions.useTimeVar=true;
    Modelssa.ssaOptions.signalUpdateRate=1;
    Modelssa.solutionScheme = 'SSA';  % Set solution scheme to SSA.
    Modelssa.ssaOptions.Nexp = 200; 
    Modelssa.ssaOptions.nSimsPerExpt = 200;
    Modelssa.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
    Modelssa.solve([],'EGRNT_SSA.csv'); 
    
    % Find MLE for simulated data set.
    EGRNT_Model_check = EGRNT_Model;  
    EGRNT_Model_check.fittingOptions.modelVarsToFit = 1:8;
    nFrePars = length(EGRNT_Model_check.fittingOptions.modelVarsToFit);
    EGRNT_MLE = zeros(1,nFrePars,Modelssa.ssaOptions.Nexp);
    fMLE = inf(1,nFrePars,Modelssa.ssaOptions.Nexp);
    B = EGRNT_Model_check;
    
    for iExp = 1:Modelssa.ssaOptions.Nexp
        
        % Dusp1  observed
        B = B.loadData('EGRNT_SSA.csv',{'x3',['exp',num2str(iExp),'_s3']}); % Link non-distorted data.
    
        fitOptions = optimset;
        fitOptions.MaxIter = 500;
        fitOptions.Display = 'iter';
        
        B.fittingOptions.timesToFit = [false,ones(1,length(B.tSpan),'logical')];
        B.tSpan = B.tSpan(2:end);
    %     B.fittingOptions.modelVarsToFit = [2,3];
         if iExp==1
             x0 = [B.parameters{B.fittingOptions.modelVarsToFit,2}]';
         else
             x0 = squeeze(EGRNT_MLE(1,:,iExp-1)); 
         end
        [EGRNT_MLE(1,:,iExp),fMLE(1,:,iExp)] = B.maximizeLikelihood(x0,fitOptions);
    end
end
    
% Solve for Fisher Information Matrix at all Time Points
B.fittingOptions.modelVarsToFit = 1:8;

B.solutionScheme = 'FSP';
B.fspOptions.fspTol = 1e-6;
B.fspOptions.bounds=[];
[fspSoln,B.fspOptions.bounds] = B.solve;

B.fspOptions.fspTol = inf;
B.solutionScheme = 'fspSens';
sensSoln_B = B.solve(fspSoln.stateSpace);

B.pdoOptions.unobservedSpecies = {'x1','x2'};
fims_B = B.computeFIM(sensSoln_B.sens);
FIM_B = B.evaluateExperiment(fims_B,B.dataSet.nCells(2:end));


% Make Plots (figure slides 14 and 15)
fimFreePars = FIM_B(B.fittingOptions.modelVarsToFit,B.fittingOptions.modelVarsToFit);
if showFigures 
    mlfeFig = figure;
    B.makeMleFimPlot(squeeze(EGRNT_MLE(1,:,:)),fimFreePars,[4,3],0.95)
    saveas(mlfeFig,'EGRNT_MLE_Fig','fig')
end

%% Using all cells to measure DUSP1 expression (smFISH)
EGRNT_Model.pdoOptions.unobservedSpecies = {'x1','x2'}; % Gene state and nuclear GR are unobserved
fims_dusp = EGRNT_Model.computeFIM(sensSoln.sens); 
FIM_dusp = EGRNT_Model.evaluateExperiment(fims_dusp,EGRNT_Model.dataSet.nCells);

nTotal_dusp = sum(EGRNT_Model.dataSet.nCells);
nCellsOpt_dusp = EGRNT_Model.optimizeCellCounts(fims_dusp,nTotal_dusp,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_dusp_case2 = EGRNT_Model.optimizeCellCounts(fims_dusp,nTotal_dusp,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt_dusp = EGRNT_Model.evaluateExperiment(fims_dusp,nCellsOpt_dusp);
fimOpt_dusp_case2 = EGRNT_Model.evaluateExperiment(fims_dusp,nCellsOpt_dusp_case2);


%% Measuring Nuclear Glucocotoroid (ICC)
EGRNT_Model.pdoOptions.unobservedSpecies = {'x3','x2'}; % Gene state and mRNA are unobserved
fims_GR = EGRNT_Model.computeFIM(sensSoln.sens);
FIM_GR = EGRNT_Model.evaluateExperiment(fims_GR,EGRNT_Model.dataSet.nCells);

nTotal_GR = sum(EGRNT_Model.dataSet.nCells);
nCellsOpt_GR = EGRNT_Model.optimizeCellCounts(fims_GR,nTotal_GR,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_GR_case2 = EGRNT_Model.optimizeCellCounts(fims_GR,nTotal_GR,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt_GR = EGRNT_Model.evaluateExperiment(fims_GR,nCellsOpt_GR);
fimOpt_GR_case2 = EGRNT_Model.evaluateExperiment(fims_GR,nCellsOpt_GR_case2);

%% Seperate experiments to measure DUSP1 and GR (separate smFISH and ICC experiments)
% The same population of cells as measuring just mRNA is used. Half of the
% cells are measuring nuclear GR while the other half measures mRNA
fims_dusp_and_gr_separate = [fims_dusp;fims_GR];

% This will fix the number of cells at half the original number for each
% experiment.
FIM = (FIM_GR+FIM_dusp)/2;

% Here we will also allow twice as many cells to be measured with this
% experiment in order to match the case above.
nTotal = sum(EGRNT_Model.dataSet.nCells);
nCellsOpt_sep = EGRNT_Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_sep_case2 = EGRNT_Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt = EGRNT_Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep);
fimOpt_case2 = EGRNT_Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep_case2);

%% Simultaneous measurement of DUSP1 and GR (smFISH and ICC experiments in the same population of cells)
EGRNT_Model.pdoOptions.unobservedSpecies = {'x2'}; % Gene state is unobserved
fims_both = EGRNT_Model.computeFIM(sensSoln.sens);
FIM_both = EGRNT_Model.evaluateExperiment(fims_both,EGRNT_Model.dataSet.nCells);

nTotal_both = sum(EGRNT_Model.dataSet.nCells);
nCellsOpt_both = EGRNT_Model.optimizeCellCounts(fims_both,nTotal_both,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_both_case2 = EGRNT_Model.optimizeCellCounts(fims_both,nTotal_both,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt_both = EGRNT_Model.evaluateExperiment(fims_both,nCellsOpt_both);
fimOpt_both_case2 = EGRNT_Model.evaluateExperiment(fims_both,nCellsOpt_both_case2);

%% Plot for Optimal Cell allocation measurement for each experiment design
timepoints = {'0min','10min','20min','30min','40min','50min','60min','75min','90min','120min','150min','180min'};

% case1 figure slide 16
figure()
designs = [EGRNT_Model.dataSet.nCells;nCellsOpt_dusp;nCellsOpt_sep(1:12);nCellsOpt_both];
designs_b = [EGRNT_Model.dataSet.nCells;nCellsOpt_dusp;nCellsOpt_sep(1:12)+nCellsOpt_sep(13:24);nCellsOpt_both];


g = bar(designs_b','barWidth',1);
hold on
h = bar(designs','barWidth',1);
hold off

set(gca,'XTickLabel',timepoints,'FontSize',20)
xlabel('Measurement timepoints')
ylabel('Number of Cell')
title('Case 1 (GR parameters known)')
legend({'','',append('Seperate GR Experiments Optimized (',int2str(sum(nCellsOpt_sep(13:24))/sum(nCellsOpt_sep)*100),'% of cells)')...
  ,'','Equal cells at each timepoint','DUSP1 Optimized',append('Seperate DUSP1 Experiments Optimized (',int2str(sum(nCellsOpt_sep(1:12))/sum(nCellsOpt_sep)*100),'% of cells)'),'Simutaneous DUSP1/GR Experiments Optimized',''},'FontSize',20)

% case2 figure slide 17
figure()
designs = [EGRNT_Model.dataSet.nCells;nCellsOpt_dusp_case2;nCellsOpt_sep_case2(1:12);nCellsOpt_both_case2];
designs_b = [EGRNT_Model.dataSet.nCells;nCellsOpt_dusp_case2;nCellsOpt_sep_case2(1:12)+nCellsOpt_sep_case2(13:24);nCellsOpt_both_case2];

g = bar(designs_b','barWidth',1);
hold on
h = bar(designs','barWidth',1);
hold off

set(gca,'XTickLabel',timepoints,'FontSize',20)
xlabel('Measurement timepoints')
ylabel('Number of Cell')
title('Case 2 (GR parameters not known)')
legend({'','',append('Seperate GR Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(13:24))/sum(nCellsOpt_sep_case2)*100),'% of cells)')...
  ,'','Equal cells at each timepoint','DUSP1 Optimized',append('Seperate DUSP1 Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(1:12))/sum(nCellsOpt_sep_case2)*100),'% of cells)'),'Simutaneous DUSP1/GR Experiments Optimized',''},'FontSize',20)

%% Comparing FIM of each experiment design (Case 1)
% figure slide 18
y = [det(inv(FIM_dusp(1:4,1:4))),det(inv(fimOpt_dusp(1:4,1:4)));...
  det(inv(FIM(1:4,1:4))),det(inv(fimOpt(1:4,1:4))); det(inv(FIM_both(1:4,1:4))), det(inv(fimOpt_both(1:4,1:4)))];
x = categorical({'Dusp1', 'Seperate GR and Dusp1 experiment', 'Simultaneous GR and Dusp1 experiment'});
bar(x,y,'barWidth',1)
set(gca, 'YScale', 'log','ylim',[1e-36,1e-15])
ax=gca;
ax.FontSize = 20;
legend({'FIM','FIM optimized'},'FontSize', 24)
ylabel('|COV|')
title('Case 1 (GR parameters known)')

%% Comparing FIM of each experiment design (Case 2)
% figure slide 19
ev = [eye(4),zeros(4)];
y = [det(ev*inv(FIM_dusp(1:8,1:8))*ev'),det(ev*inv(fimOpt_dusp_case2(1:8,1:8))*ev');...
  det(ev*inv(FIM(1:8,1:8))*ev'),det(ev*inv(fimOpt_case2(1:8,1:8))*ev'); det(ev*inv(FIM_both(1:8,1:8))*ev'), det(ev*inv(fimOpt_both_case2(1:8,1:8))*ev')];
x = categorical({'Dusp1', 'Seperate GR and Dusp1 experiment', 'Simultaneous GR and Dusp1 experiment'});
bar(x,y,'barWidth',1)
set(gca, 'YScale', 'log','ylim',[1e-25,1e-8])
ax2=gca;
ax2.FontSize = 25;
legend({'FIM','FIM optimized'},'FontSize', 24)
ylabel('|COV|')
title('Case 2 (GR parameters not known')
nTotal_all = sum(EGRNT_Model.dataSet.nCells);

%% Effect of additional cell measurements on information gained
% figure slide 21
y = [det(inv(FIM(1:4,1:4))),det(inv(fimOpt_dusp(1:4,1:4)))];
x = categorical({'Simple DUSP1 Model'});
if showFigures
    figure()
    bar(x,y,'barWidth',1)
    set(gca, 'YScale', 'log','ylim',[1e-36,1e-15])
    ax=gca;
    ax.FontSize = 20;
    legend({'FIM','FIM optimized'},'FontSize', 24)
    ylabel('|FIM^-^1|')
    title('Simple DUSP1 Model Optimization')
end
%% Tryptolide Experiment
if runOptTriptolideExp

    ModelTrypt = EGRNT_Model;
    
    % Changing transcription propensity function to include transcription repression mechanism
    ModelTript.propensityFunctions(5) = {'kr*x2*Itrypt'};
    ModelTript.inputExpressions(2,:) = {'Itrypt','(t<tpt)'}; % Repression dynamic occurs after time equal to 'tpt'
    
    ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';
    
    tpt_array = [0:5:40,42:2:180]; % Range of times to calculate FIM  
    sensSoln = cell(1,length(tpt_array));
    
    % Solve the FIM at each time point of the experiment to verify the optimal
    % time to apply triptolide in an experiment
    for iTpt = 1:length(tpt_array)
    
        ModelTrypt.parameters(8,:) = {'tpt',tpt_array(iTpt)};
        
        % Get FSP fit for bounds.
        ModelTrypt.solutionScheme = 'FSP';
        ModelTrypt.fspOptions.fspTol = 1e-6;
        ModelTrypt.fspOptions.bounds=[];
        [fspSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;
        
        % Calculate Sensitivity Matrix from FSP solution
        ModelTrypt.fspOptions.fspTol = inf;
        ModelTrypt.solutionScheme = 'fspSens';
        sensSoln{iTpt} = ModelTrypt.solve(fspSoln.stateSpace);
        
        % Solve for Fisher Information Matrix at all Time Points using sensitivity matrix
        ModelTrypt.pdoOptions.unobservedSpecies = {'x1','x2'};
        fims = ModelTrypt.computeFIM(sensSoln{iTpt}.sens);
        FIM = ModelTrypt.evaluateExperiment(fims,ModelTrypt.dataSet.nCells);
        expectedDetCov(iTpt) = det(FIM^(-1));
    
    end
    
    % |FIM^-1| vs triptilde application time
    % figure slide 25
    if figure
        figure;plot(tpt_array,expectedDetCov); hold on
        set(gca,'yscale','log')
        ylim = get(gca,'ylim');
        for it = 1:length(ModelTrypt.dataSet.times)
        plot(ModelTrypt.dataSet.times(it)*[1,1],ylim,'k--')
        end
        legend('|FIM^-1|','Original Experiment Measurement Times')
        xlabel('Time [min]')
        ylabel('|FIM^-1|')
    end

end

%% Save Data
if saveMHOutput % Save mhResults
    i=1;
        while exist(['EGRNT_mhast','_v',num2str(i),'.mat'],'file')
            i=i+1;
        end
    save(['EGRNT_mhast','_v',num2str(i),'.mat'],'mhResults');
end

if saveModelOutput % Save Model
    j=1;
        while exist(['EGRNT_Model','_v',num2str(j),'.mat'],'file')
            j=j+1;
        end
    save(['EGRNT_Model','_v',num2str(j),'.mat'],'EGRNT_Model');

end

if saveSensOutput   % Save Sensitivity Results
    k=1;
        while exist(['EGRNT_sens','_v',num2str(k),'.mat'],'file')
            k=k+1;
        end
    save(['EGRNT_sens','_v',num2str(k),'.mat'],'sensSoln');

end

if saveMLEOutput   % Save MLE Results
    n=1;
        while exist(['EGRNT_MLE','_v',num2str(n),'.mat'],'file')
            n=n+1;
        end
    save(['EGRNT_MLE','_v',num2str(n),'.mat'],'EGRNT_MLE');
    save(['B.mat'],'B') % Save model that generated the MLE

end

if saveTriptOutput   % Save potimal Tript experiment Results
    n=1;
        while exist(['EGRNT_Tript_FIM','_v',num2str(n),'.mat'],'file')
            n=n+1;
        end
    save(['EGRNT_Tript_FIM','_v',num2str(n),'.mat'],'FIM_Tript');

end

