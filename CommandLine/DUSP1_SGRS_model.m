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

calcMLE = true; % MLE calculation for FIM verification
runOptTriptolideExp = false; % optimal triptolide application calc

showFigures = false;

% Save Model, MHresults, sens, and MLE results
% Previous results will not be overwritten
saveMHOutput = true; 
saveModelOutput = true; 
saveSensOutput = true; 
saveMLEOutput = true; 
saveTriptOutput = false;

%% Define SSIT Model
% Initial Set up for the model
Model = SSIT; % Rename SSIT code to make changes with out changing the original code
Model.species = {'x1';'x2'}; % x1: number of on alleles,  x2: mRNA 
Model.initialCondition = [0;0]; % Set initial conditions
Model.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'}; % Set the propensity functions of the 4 reactions

% Associated stoichiometry of the four reactions
stoich = [1,0;   % Gene allele switching to on state
         -1,0;   % Gene allele switching to off state
          0,1;   % mRNA production
          0,-1]; % mRNA degredation

Model.stoichiometry = transpose(stoich);
% Input expression for time varying signals
Model.inputExpressions = {'IGR','kcn0/knc+(t>=0)*kcn1/(r1-knc)*(exp(-knc*t)-exp(-r1*t))'};

% Defining the values of each parameter used
Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
                   'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});

Model.fspOptions.initApproxSS = true;


SGRS_Model = Model;
tpt_array = 0:20:180;
SGRS_Model.tSpan = tpt_array; % Define the time points

%% Load saved mHast and Sensitivity
if useSavedModel
    j=1;
    while exist(['SGRS_model','_v',num2str(j),'.mat'],'file')
        j=j+1;
    end

    if exist(['SGRS_model','_v',num2str(j-1),'.mat'],'file')
        SGRS_Model = load('simple_dusp1_model.mat').simple_Model;
%         Model = load('simple_dusp1_model.mat').simple_Model;
    else
        error(['No file named ','SGRS_model','_v',num2str(j-1),'.mat found'])
    end

end

if useSavedModelMH
    i=1;
    while exist(['SGRS_mhast','_v',num2str(i),'.mat'],'file')
        i=i+1;
    end
    if exist(['SGRS_mhast','_v',num2str(i-1),'.mat'],'file')
        SGRS_Model = load(['SGRS_mhast','_v',num2str(i-1),'.mat']).mhResults;
        %     mhResults = load('complex_dusp1_mhast.mat').mhResults;
    else
        error(['No file named ','SGRS_mhast','_v',num2str(i),'.mat',' file found'])
    end

end

if useSavedModelSensitivty
    k=1;
    while exist(['SGRS_sens','_v',num2str(k),'.mat'],'file')
        k=k+1;
    end
    if exist(['SGRS_sens','_v',num2str(k-1),'.mat'],'file')
        sensSoln = load(['SGRS_sens','_v',num2str(k-1),'.mat']).sensSoln;
%     sensSoln = load('complex_dusp1_sens.mat').sensSoln;
    else
        error(['No file named ','SGRS_sens','_v',num2str(k-1),'.mat found'])
    end
end

if useSavedModelMLE
    n=1;
    while exist(['SGRS_MLE','_v',num2str(n),'.mat'],'file')
        n=n+1;
    end
    if exist(['SGRS_MLE','_v',num2str(n-1),'.mat'],'file')
        SGRS_MLE = load(['SGRS_MLE','_v',num2str(n-1),'.mat']).sensSoln;
    else
        error(['No file named ','SGRS_MLE','_v',num2str(n),'.mat found'])
    end
end
%% Solve the model using the FSP
SGRS_Model = SGRS_Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'}); % Load experimental data set
if fitModel
    SGRS_Model.fittingOptions.modelVarsToFit = 1:8; % Picking parameters 1-8 for fitting
    
    for i=1:5
        SGRS_Model.solutionScheme = 'FSP';
        SGRS_Model.fspOptions.fspTol = 1e-4;
        SGRS_Model.fspOptions.verbose = 0;
        SGRS_Model.fspOptions.bounds=[];
        [fspSoln,SGRS_Model.fspOptions.bounds] = SGRS_Model.solve;
        SGRS_Model.fspOptions.bounds
        
        % Load and Fit smFISH Data
        SGRS_Model.fspOptions.fspTol = inf;
        fitOptions = optimset('Display','iter','MaxIter',500); 
        SGRS_Model.parameters(SGRS_Model.fittingOptions.modelVarsToFit,2) =...
          num2cell(SGRS_Model.maximizeLikelihood(...
          [SGRS_Model.parameters{SGRS_Model.fittingOptions.modelVarsToFit,2}],...
          fitOptions));
    end
    
    % figure slide 7 and 8
    if showFigures
        SGRS_Model.makeFitPlot;
    end
end
%% Metropolis Hastings to Quantify Parameter Uncertainty
if runMetHast
    SGRS_Model.fittingOptions.modelVarsToFit = 1:4; % Picking parameters 1-4 for fitting
    
    allFitOptions.CovFIMscale = 0.1;% make smaller for higher acceptance rate
    MHOptions = struct('numberOfSamples',15000,'burnin',100,'thin',1,...
      'useFIMforMetHast',true,'suppressFSPExpansion',true);
    [bestParsFound,~,mhResults] = SGRS_Model.maximizeLikelihood([SGRS_Model.parameters{SGRS_Model.fittingOptions.modelVarsToFit,2}]',...
      MHOptions,'MetropolisHastings');
    SGRS_Model.parameters(SGRS_Model.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
    
    % figure slide 10
    if showFigures
        SGRS_Model.plotMHResults(mhResults); 
    end
end

%% Calculate CME Sensitivity to Parameter Variations
SGRS_Model.sensOptions.solutionMethod = 'finiteDifference';
SGRS_Model.solutionScheme = 'fspSens'; % Choosing sensitivity solution method
SGRS_Model.fspOptions.fspTol = 1e-6;
[sensSoln,SGRS_Model.fspOptions.bounds] = SGRS_Model.solve;

%% Solve for Fisher Information Matrix at all Time Points
SGRS_Model.pdoOptions.unobservedSpecies = {'x1'}; % PDO applied for the case where the Gene state is not observed

% FIM calculation taking into account the number of cells measured at each time point
fims = SGRS_Model.computeFIM(sensSoln.sens);
FIM = SGRS_Model.evaluateExperiment(fims,SGRS_Model.dataSet.nCells);

nTotal_dusp = sum(SGRS_Model.dataSet.nCells);
nCellsOpt_dusp = SGRS_Model.optimizeCellCounts(fims,nTotal_dusp,'TR[1:4]');
fimOpt_dusp = SGRS_Model.evaluateExperiment(fims,nCellsOpt_dusp); % Using the FIM to optimize the cell measuremet times

if showFigures
    SGRS_Model.plotMHResults(mhResults,fimOpt_dusp);
end

%% Maximum likelihood Estimator Verification
if calcMLE
    % Simmulate data with known parameters and re-fit simmulated data to get MLE 
    starttime=-1500;
    SGRSssa=SGRS_Model;
    SGRSssa.initialTime = starttime;
    SGRSssa.tSpan = [starttime,SGRS_Model.tSpan];
    SGRSssa.ssaOptions.useTimeVar=true;
    SGRSssa.ssaOptions.signalUpdateRate=1;
    SGRSssa.solutionScheme = 'SSA';  % Set solution scheme to SSA.
    SGRSssa.ssaOptions.Nexp = 200; 
    SGRSssa.ssaOptions.nSimsPerExpt = 200;
    SGRSssa.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
    SGRSssa.solve([],'SGRS_SSA.csv'); 
    
    % Find MLE for simulated data set.
    
    SGRS_Model_check.fittingOptions.modelVarsToFit = 1:8;
    nFrePars = length(SGRS_Model_check.fittingOptions.modelVarsToFit);
    SGRS_MLE = zeros(1,nFrePars,SGRSssa.ssaOptions.Nexp);
    fMLE = inf(1,nFrePars,SGRSssa.ssaOptions.Nexp);
    B = SGRS_Model_check.loadData('SGRS_SSA.csv',{'x1','exp1_s1';'x2','exp1_s2'});
    
    for iExp = 1:SGRSssa.ssaOptions.Nexp
        
        % Dusp1  observed
        B = SGRS_Model_check.loadData('SGRS_SSA.csv',{'x2',['exp',num2str(iExp),'_s2']}); % Link non-distorted data.
    
        fitOptions = optimset;
        fitOptions.MaxIter = 500;
        fitOptions.Display = 'iter';
        
        B.fittingOptions.timesToFit = [false,ones(1,length(SGRS_Model_check.tSpan),'logical')];
        B.tSpan = B.tSpan(2:end);
    %     B.fittingOptions.modelVarsToFit = [2,3];
         if iExp==1
             x0 = [B.parameters{B.fittingOptions.modelVarsToFit,2}]';
         else
             x0 = squeeze(SGRS_MLE(1,:,iExp-1)); 
         end
        [SGRS_MLE(1,:,iExp),fMLE(1,:,iExp)] = B.maximizeLikelihood(x0,fitOptions);
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
    
    B.pdoOptions.unobservedSpecies = {'x1'};
    fims_B = B.computeFIM(sensSoln_B.sens);
    FIM_B = B.evaluateExperiment(fims_B,B.dataSet.nCells);
    
    
    % Make Plots (figure slides 14 and 15)
    fimFreePars = FIM_B(B.fittingOptions.modelVarsToFit,B.fittingOptions.modelVarsToFit);
    if showFigures    
        B.makeMleFimPlot(squeeze(SGRS_MLE(1,:,:)),fimFreePars,[4,3],0.95)
    end
end

%%  Plot cell's per time stack
timepoints = {'0min','10min','20min','30min','40min','50min','60min','75min','90min','120min','150min','180min'};

% figure slide 20
if showFigures
    figure()
    designs = [SGRS_Model.dataSet.nCells';nCellsOpt_dusp;];
    
    g = bar(designs','barWidth',1);
    
    set(gca,'XTickLabel',timepoints,'FontSize',20)
    xlabel('Measurement timepoints')
    ylabel('Number of Cell')
    title('Simple DUSP1 Model Optimization')
    legend({'Equal Number of cells at each timepoint','Optimized experiment'});
end
%% FIM bar plots 
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

    ModelTrypt = SGRS_Model;
    
    % Changing transcription propensity function to include transcription repression mechanism
    ModelTrypt.propensityFunctions(3) = {'kr*x1*Itrypt'};
    ModelTrypt.inputExpressions(2,:) = {'Itrypt','(t<tpt)'}; % Repression dynamic occurs after time equal to 'tpt'
    
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
        ModelTrypt.pdoOptions.unobservedSpecies = {'x1'};
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
        while exist(['SGRS_mhast','_v',num2str(i),'.mat'],'file')
            i=i+1;
        end
    save(['SGRS_mhast','_v',num2str(i),'.mat'],'mhResults');
end

if saveMHOutput % Save Model
    j=1;
        while exist(['SGRS_model','_v',num2str(j),'.mat'],'file')
            j=j+1;
        end
    save(['SGRS_model','_v',num2str(j),'.mat'],'SGRS_Model');

end

if saveSensOutput   % Save Sensitivity Results
    k=1;
        while exist(['SGRS_sens','_v',num2str(k),'.mat'],'file')
            k=k+1;
        end
    save(['SGRS_sens','_v',num2str(k),'.mat'],'sensSoln');

end

if saveMLEOutput   % Save MLE Results
    n=1;
        while exist(['SGRS_MLE','_v',num2str(n),'.mat'],'file')
            n=n+1;
        end
    save(['SGRS_MLE','_v',num2str(n),'.mat'],'SGRS_MLE');

end

if saveTriptOutput   % Save potimal Tript experiment Results
    n=1;
        while exist(['SGRS_Tript_FIM','_v',num2str(n),'.mat'],'file')
            n=n+1;
        end
    save(['SGRS_Tript_FIM','_v',num2str(n),'.mat'],'FIM_Tript');

end

