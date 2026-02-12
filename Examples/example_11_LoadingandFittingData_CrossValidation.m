%% SSIT/Examples/example_11_LoadingandFittingData_CrossValidation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Fitted model evaluation using cross-validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the 4-state STL1 model from example_1_CreateSSITModels and its FSP 
% solutions from example_4_SolveSSITModels_FSP
%clear
%close all
addpath(genpath('../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE
% example_10_LoadingandFittingData_MH

%% Load model fitted using Metropolis-Hastings:
% load('example_10_LoadingandFittingData_MH.mat')

% View summary of 4-state STL1 model:
STL1_4state_MH.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Cross-Validation fitting different replicas at the same time
%  using SSITMultiModel
%
% In this example, we will use the multimodel to allow parameters to change
% for different replica data sets (e.g., to allow for batch variations, or
% to explore how parameters change under different genetic variations).
% Here, we illustrate a quick means to generate the multimodel starting
% with a single template and a datafile with multiple replicas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Fitting Options:
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display','final','MaxIter',200);

% Make a copy of our 4-state STL1 model:
STL1_4state_CrossVal = STL1_4state_MH;

% Specify datafile name and species linking rules:
DataFileName = 'data/filtered_data_2M_NaCl_Step.csv';
LinkedSpecies = {'mRNA','RNA_STL1_total_TS3Full'};

% Suppose we only wish to fit the data at times before 25 minutes.  
% Set the global conditions:
ConditionsGlobal = {[],[],'TAB.time<=25'};

% Split up the replicas to be separate:
ConditionsReplicas = {'TAB.Replica==1';'TAB.Replica==2'};

% Specify constraints on rep-to-rep parameter variations. Here, we specify 
% that there is an expected 0.1 log10 deviation expected in some parameters 
% and smaller in others.  No deviation at all is indicated by 0.
Log10Constraints = ...
    [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.02,0.02,0.02,0.02,0.02]; 

% Create full model:
CrossValidationModel = SSITMultiModel.createCrossValMultiModel(...
    STL1_4state_CrossVal, DataFileName, LinkedSpecies, ConditionsGlobal,...
    ConditionsReplicas, Log10Constraints);
CrossValidationModel = CrossValidationModel.initializeStateSpaces;

% Run the model fitting routines:
crossValPars = CrossValidationModel.parameters;
crossValPars = CrossValidationModel.maximizeLikelihood(...
    crossValPars, fitOptions, fitAlgorithm);
CrossValidationModel = CrossValidationModel.updateModels(crossValPars);
CrossValidationModel.parameters = crossValPars;

% Get the fixed parameters (Hog1 input):
fixed = cell2mat(STL1_4state_CrossVal.parameters(14:18,2)).';  % 1x5 row

% Combine newly fit parameters from each replica model with fixed 
% parameters (1x36 vector):
CrossValidationModel.parameters = ...
    [CrossValidationModel.parameters(1:13), fixed,...
     CrossValidationModel.parameters(14:26), fixed];    

% Update parameter indices:
CrossValidationModel.parameterIndices = {1:18, 19:36};

% Make a figure to explore how much the parameters changed between replicas:
fignum = 12; useRelative = true;
CrossValidationModel.compareParameters(fignum,useRelative);

%% Save model & cross-validation results:
saveNames = unique({ ...
    'STL1_4state_CrossVal'
    'ConditionsGlobal'
    'ConditionsReplicas'
    'crossValPars'
    'CrossValidationModel'
    });
    
save('example_11_LoadingandFittingData_CrossValidation',saveNames{:})