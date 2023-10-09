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
saveMLEOutput = false; 
saveTriptOutput = false;

%% Define SSIT Model (Glucocortoroid Receptor Translocation)
% Initial Set up for the model
ModelGR = SSIT;
ModelGR.species = {'cytGR';'nucGR'};
ModelGR.initialCondition = [20;1];

% Associated stoichiometry of the four reactions
ModelGR.propensityFunctions = {'(kcn0+kcn1*IDex)*cytGR';'knc*nucGR';'kg1';'gg1*cytGR'};

stoich = [-1,1; % GR translocates to nucleus
          1,-1; % GR translocates out of nucleus
          1,0;  % Cytoplasm GR creation
          -1,0];% Cytoplasm GR degredation
ModelGR.stoichiometry = transpose(stoich);

ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
    'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',0.0000;'gg1',0.00000;'gg2',0.0001});

Model.fspOptions.initApproxSS = true;

ModelGR.dataSet = [];

GR_Model = Model;
GR_Model.inputExpressions = {'IDex','100*exp(-gDex*t)*(t>0)'}; % Model for 100nM Dex stimulation

tpt_array = 0:20:180;
GR_Model.tSpan = tpt_array; % Define the time points
ModelGR.fittingOptions.modelVarsToFit = [5,6,7,8,9,10,11];

%% Load saved mHast and Sensitivity
if useSavedModel
    j=1;
    while exist(['GR_model','_v',num2str(j),'.mat'],'file')
        j=j+1;
    end

    if exist(['GR_model','_v',num2str(j-1),'.mat'],'file')
        GR_Model = load(['GR_model','_v',num2str(j-1),'.mat']).GR_Model;
%         Model = load('simple_dusp1_model.mat').simple_Model;
    else
        error(['No file named ','GR_model','_v',num2str(j-1),'.mat found'])
    end

end

if useSavedModelMH
    i=1;
    while exist(['GR_mhast','_v',num2str(i),'.mat'],'file')
        i=i+1;
    end
    if exist(['GR_mhast','_v',num2str(i-1),'.mat'],'file')
        mhResults = load(['GR_mhast','_v',num2str(i-1),'.mat']).mhResults;
        %     mhResults = load('complex_dusp1_mhast.mat').mhResults;
    else
        error(['No file named ','GR_mhast','_v',num2str(i),'.mat',' file found'])
    end

end

if useSavedModelSensitivty
    k=1;
    while exist(['GR_sens','_v',num2str(k),'.mat'],'file')
        k=k+1;
    end
    if exist(['GR_sens','_v',num2str(k-1),'.mat'],'file')
        sensSoln = load(['GR_sens','_v',num2str(k-1),'.mat']).sensSoln;
%     sensSoln = load('complex_dusp1_sens.mat').sensSoln;
    else
        error(['No file named ','GR_sens','_v',num2str(k-1),'.mat found'])
    end
end

if useSavedModelMLE
    n=1;
    while exist(['GR_MLE','_v',num2str(n),'.mat'],'file')
        n=n+1;
    end
    if exist(['GR_MLE','_v',num2str(n-1),'.mat'],'file')
        GR_MLE = load(['GR_MLE','_v',num2str(n-1),'.mat']).GR_MLE;
        B = load(['B.mat']).B;
    else
        error(['No file named ','GR_MLE','_v',num2str(n),'.mat found'])
    end
end
%% Solve the model using the FSP
ModelGR100 = ModelGR100.loadData('DUSP1_GR_dataframes/GR_ICC_3hr_Dex_total.csv',...
    {'nucGR','GRNorm'},{'conc','100'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

if fitModel
    
    for i=1:5
        GR_Model.solutionScheme = 'FSP';
        GR_Model.fspOptions.fspTol = 1e-4;
        GR_Model.fspOptions.verbose = 0;
        GR_Model.fspOptions.bounds=[];
        [fspSoln,GR_Model.fspOptions.bounds] = GR_Model.solve;
        GR_Model.fspOptions.bounds
        
        % Load and Fit smFISH Data
        GR_Model.fspOptions.fspTol = inf;
        fitOptions = optimset('Display','iter','MaxIter',500); 
        GR_Model.parameters(GR_Model.fittingOptions.modelVarsToFit,2) =...
          num2cell(GR_Model.maximizeLikelihood(...
          [GR_Model.parameters{GR_Model.fittingOptions.modelVarsToFit,2}],...
          fitOptions));
    end
    
    % figure slide 7 and 8
    if showFigures
        GR_Model.makeFitPlot;
    end
end
%% Metropolis Hastings to Quantify Parameter Uncertainty
if runMetHast
    
    allFitOptions.CovFIMscale = 0.1;% make smaller for higher acceptance rate
    MHOptions = struct('numberOfSamples',15000,'burnin',100,'thin',1,...
      'useFIMforMetHast',true,'suppressFSPExpansion',true);
    [bestParsFound,~,mhResults] = GR_Model.maximizeLikelihood([GR_Model.parameters{GR_Model.fittingOptions.modelVarsToFit,2}]',...
      MHOptions,'MetropolisHastings');
    GR_Model.parameters(GR_Model.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
    
    % figure slide 10
    if showFigures
        GR_Model.plotMHResults(mhResults); 
    end
end

%% Calculate CME Sensitivity to Parameter Variations
GR_Model.sensOptions.solutionMethod = 'finiteDifference';
GR_Model.solutionScheme = 'fspSens'; % Choosing sensitivity solution method
GR_Model.fspOptions.fspTol = 1e-6;
[sensSoln,GR_Model.fspOptions.bounds] = GR_Model.solve;

%% Solve for Fisher Information Matrix at all Time Points
GR_Model.pdoOptions.unobservedSpecies = {'x1'}; % PDO applied for the case where the Gene state is not observed

% FIM calculation taking into account the number of cells measured at each time point
fims = GR_Model.computeFIM(sensSoln.sens);
FIM = GR_Model.evaluateExperiment(fims,GR_Model.dataSet.nCells);

nTotal_dusp = sum(GR_Model.dataSet.nCells);
nCellsOpt_dusp = GR_Model.optimizeCellCounts(fims,nTotal_dusp,'TR[1:4]');
fimOpt_dusp = GR_Model.evaluateExperiment(fims,nCellsOpt_dusp); % Using the FIM to optimize the cell measuremet times

if showFigures
    GR_Model.plotMHResults(mhResults,fimOpt_dusp);
end

%% Maximum likelihood Estimator Verification
if calcMLE
    % Simmulate data with known parameters and re-fit simmulated data to get MLE 
    starttime=-1500;
    GRssa=GR_Model;
    GRssa.initialTime = starttime;
    GRssa.tSpan = [starttime,GR_Model.tSpan];
    GRssa.ssaOptions.useTimeVar=true;
    GRssa.ssaOptions.signalUpdateRate=1;
    GRssa.solutionScheme = 'SSA';  % Set solution scheme to SSA.
    GRssa.ssaOptions.Nexp = 200; 
    GRssa.ssaOptions.nSimsPerExpt = 200;
    GRssa.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
    GRssa.solve([],'GR_SSA.csv'); 
    
    % Find MLE for simulated data set.
    GR_Model_check = GR_Model;  
    GR_Model_check.fittingOptions.modelVarsToFit = [5,6,7,8,9,10,11];
    nFrePars = length(GR_Model_check.fittingOptions.modelVarsToFit);
    GR_MLE = zeros(1,nFrePars,GRssa.ssaOptions.Nexp);
    fMLE = inf(1,nFrePars,GRssa.ssaOptions.Nexp);
%     B = GR_Model_check.loadData('GR_SSA.csv',{'x1','exp1_s1';'x2','exp1_s2'});
    B = GR_Model_check;
    
    for iExp = 1:GRssa.ssaOptions.Nexp
        
        % Dusp1  observed
        B = B.loadData('GR_SSA.csv',{'x2',['exp',num2str(iExp),'_s2']}); % Link non-distorted data.
    
        fitOptions = optimset;
        fitOptions.MaxIter = 500;
        fitOptions.Display = 'iter';
        
        B.fittingOptions.timesToFit = [false,ones(1,length(B.tSpan),'logical')];
        B.tSpan = B.tSpan(2:end);

        if iExp==1
             x0 = [B.parameters{B.fittingOptions.modelVarsToFit,2}]';
         else
             x0 = squeeze(GR_MLE(1,:,iExp-1)); 
         end
        [GR_MLE(1,:,iExp),fMLE(1,:,iExp)] = B.maximizeLikelihood(x0,fitOptions);
    end
end
    
% Solve for Fisher Information Matrix at all Time Points
B.fittingOptions.modelVarsToFit = [5,6,7,8,9,10,11];

B.solutionScheme = 'FSP';
B.fspOptions.fspTol = 1e-6;
B.fspOptions.bounds=[];
[fspSoln,B.fspOptions.bounds] = B.solve;

B.fspOptions.fspTol = inf;
B.solutionScheme = 'fspSens';
sensSoln_B = B.solve(fspSoln.stateSpace);

B.pdoOptions.unobservedSpecies = {'x1'};
fims_B = B.computeFIM(sensSoln_B.sens);
FIM_B = B.evaluateExperiment(fims_B,B.dataSet.nCells(2:end));


% Make Plots (figure slides 14 and 15)
fimFreePars = FIM_B(B.fittingOptions.modelVarsToFit,B.fittingOptions.modelVarsToFit);
if showFigures 
    mlfeFig = figure;
    B.makeMleFimPlot(squeeze(GR_MLE(1,:,:)),fimFreePars,[4,3],0.95)
    saveas(mlfeFig,'GR_MLE_Fig','fig')
end



if runOptTriptolideExp

    ModelTrypt = GR_Model;
    
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
        while exist(['GR_mhast','_v',num2str(i),'.mat'],'file')
            i=i+1;
        end
    save(['GR_mhast','_v',num2str(i),'.mat'],'mhResults');
end

if saveModelOutput % Save Model
    j=1;
        while exist(['GR_model','_v',num2str(j),'.mat'],'file')
            j=j+1;
        end
    save(['GR_model','_v',num2str(j),'.mat'],'GR_Model');

end

if saveSensOutput   % Save Sensitivity Results
    k=1;
        while exist(['GR_sens','_v',num2str(k),'.mat'],'file')
            k=k+1;
        end
    save(['GR_sens','_v',num2str(k),'.mat'],'sensSoln');

end

if saveMLEOutput   % Save MLE Results
    n=1;
        while exist(['GR_MLE','_v',num2str(n),'.mat'],'file')
            n=n+1;
        end
    save(['GR_MLE','_v',num2str(n),'.mat'],'GR_MLE');
    save(['B.mat'],'B') % Save model that generated the MLE

end

