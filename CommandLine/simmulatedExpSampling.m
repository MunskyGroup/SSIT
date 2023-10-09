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
%% Run simulated experiment sampling times
timeMatrix = [0,90,180];
NCells = 100;
dataFileName =  'DUSP1_3hr_Dex_100nM_total.csv'; 
[simData,csvFile] = sampleExperimentSim(dataFileName,timeMatrix,NCells);

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

tpt_array = [0 10 20 30 40 50 60 75 90 120 150 180];
Model.tSpan = tpt_array; % Define the time points

%% Solve the model using the FSP
Model.fittingOptions.modelVarsToFit = 1:8; % Picking parameters 1-8 for fitting

% Priors for fitting
muLog10Prior = [-1 -1 1 -2 -1 -2 -1 -1];
sigLog10Prior = [2 2 1 1 2 2 2,2];

muLog10Prior = muLog10Prior(Model.fittingOptions.modelVarsToFit);
sigLog10Prior = sigLog10Prior(Model.fittingOptions.modelVarsToFit);
Model.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);

Model = Model.loadData(csvFile,{'x2','RNA_nuc'}); % Load experimental data set

for i=1
    Model.solutionScheme = 'FSP';
    Model.fspOptions.fspTol = 1e-4;
    Model.fspOptions.verbose = 0;
    Model.fspOptions.bounds=[];
    [fspSoln,Model.fspOptions.bounds] = Model.solve;
    Model.fspOptions.bounds
    
    % Load and Fit smFISH Data
    Model.fspOptions.fspTol = inf;
    fitOptions = optimset('Display','iter','MaxIter',500); 
    Model.parameters(Model.fittingOptions.modelVarsToFit,2) =...
      num2cell(Model.maximizeLikelihood(...
      [Model.parameters{Model.fittingOptions.modelVarsToFit,2}],...
      fitOptions));
end
Model.tSpan
Model.makeFitPlot
%% Metropolis Hastings to Quantify Parameter Uncertainty
Model.fittingOptions.modelVarsToFit = 1:4; % Picking parameters 1-4 for fitting
muLog10Prior = [-1 -1 1 -2];
sigLog10Prior = [2 2 1 1];

muLog10Prior = muLog10Prior(Model.fittingOptions.modelVarsToFit);
sigLog10Prior = sigLog10Prior(Model.fittingOptions.modelVarsToFit);
Model.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
Model.fittingOptions.logPrior = []; % Remove prior before fitting

allFitOptions.CovFIMscale = 0.1;% make smaller for higher acceptance rate
MHOptions = struct('numberOfSamples',5000,'burnin',100,'thin',1,...
  'useFIMforMetHast',true,'suppressFSPExpansion',true);
[bestParsFound,~,mhResults] = Model.maximizeLikelihood([Model.parameters{Model.fittingOptions.modelVarsToFit,2}]',...
  MHOptions,'MetropolisHastings');
Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
Model.tSpan
%% Pick Samples
mhSamplePar = exp(mhResults.mhSamples);
pickN = 20; % Number of samples to pick
halfStack = length(mhSamplePar(:,1))/2;
k1 = randperm(halfStack,pickN);
k2 = halfStack + k1; % Picks random samples from second half of stack
par_pick = mhSamplePar(k2,:);
fims_sim = zeros(pickN,length(tpt_array,Model.parameters),length(Model.parameters));

for i=1:pickN
    ModelSamplePar = Model;
    ModelSamplePar.parameters = par_pick(i,:);
    ModelSamplePar.parameters = ({'koff',par_pick(i,1);'kon',par_pick(i,2);'kr',par_pick(i,3);'gr',par_pick(i,4);...
                   'kcn0',Model.parameters{5,2};'kcn1',Model.parameters{6,2};'knc',Model.parameters{7,2};'r1',Model.parameters{8,2}});

%     ModelSamplePar.tSpan = [0,3,6,...]
    ModelSamplePar.sensOptions.solutionMethod = 'finiteDifference';
    ModelSamplePar.solutionScheme = 'fspSens'; % Choosing sensitivity solution method
    ModelSamplePar.fspOptions.fspTol = 1e-6;
    [sensSoln,ModelSamplePar.fspOptions.bounds] = ModelSamplePar.solve;

    SGRS_Model.pdoOptions.unobservedSpecies = {'x1'}; % PDO applied for the case where the Gene state is not observed
    
    % FIM calculation taking into account the number of cells measured at each time point
    fims_sim(i,:,:,:) = Model.computeFIM(sensSoln.sens);
%     FIM_sim(i,:,:) = Model.evaluateExperiment(fims_sim(i,:,:,:),Model.dataSet.nCells);
end

fim_average = squeeze(mean(fims_sim,1));



%% Average FIM of Samples
FIM_Sum = zeros(length(Model.parameters(:,1)),length(Model.parameters(:,1)));
for i = 1:length(FIM_sim(:,:,1))
    FIM_Sum = FIM_Sum + squeeze(FIM_sim(i,:,:))^(-1);
end
FIM_inv_avg = FIM_Sum/length(FIM_sim(:,:,1));
COV = det(FIM_inv_avg);

%% Suggested New Measurement time
% addNewTime = [10:10:180];
% for i = 1:length(addNewTime)
%     newTimeMatrix = timeMatrix;
%     newTimeMatrix = [newTimeMatrix, addNewTime(i)];
%     newTimeMatrix = sort(newTimeMatrix);
% 
%     dataFileName =  'DUSP1_3hr_Dex_100nM_total.csv'; 
%     [simData,csvFile2] = sampleExperimentSim(dataFileName,newTimeMatrix,NCells);
%     Model.tSpan = newTimeMatrix;
%     Model.loadData(csvFile2,{'x2','RNA_nuc'})
% end