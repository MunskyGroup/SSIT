%% example_11_LoadingandFittingData_CrossValidation

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

%% Load pre-computed FSP solutions:
% load('example_4_SolveSSITModels_FSP.mat')

% View summary of 4-state STL1 model:
STL1_4state_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform Cross-Validation fitting different replicas at the same time
%  using SSITMultiModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of our 4-state STL1 model for Metropolis-Hastings (MH):
STL1_4state_CrossVal = STL1_4state_FSP;

% Specify datafile name and species linking rules
DataFileName = 'data/filtered_data_2M_NaCl_Step.csv';
LinkedSpecies = {'mRNA','RNA_STL1_total_TS3Full'};

% In this case, let's suppose that we only wish to fit the data at times
% before 75 minutes.  We will set the global conditions as:
ConditionsGlobal = {[],[],'TAB.time<=75'};

% We want to split up the replicas to be separate.
ConditionsReplicas = {'TAB.Replica==1';...
    'TAB.Replica==2'};

modelLibrary = 'crossValidationModels_STL1';

STL1_4state_CrossVal.runCrossValidation(STL1_4state_CrossVal,...
    DataFileName,LinkedSpecies,ConditionsGlobal,ConditionsReplicas, ...
    modelLibrary)

