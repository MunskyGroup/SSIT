%% example_11_ComplexModels_MultiModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Fitting multiple models and data sets with shared parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels 

% View model summaries:
STL1_data.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit Multiple Models and Data sets with Shared Parameters
%  Example script to show how multiple SSIT models and data sets can be 
%  fit simultaneously.  This is most useful in situations where:
%     1) the analysis considers different experimental conditions (e.g.,
%        different time points, different inducer concentrations, 
%        different genetic mutations);
%     2) replica to replica variations are expected that would result in
%        slightly different parameter combinations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the STL1 model for Replica 1:
STL1_data_rep1 = STL1_data;

% Set FSP tolerance and model parameters to fit:
STL1_data_rep1.fspOptions.fspTol = inf;
STL1_data_rep1.fittingOptions.modelVarsToFit = 1:7;

%% Load and associate smFISH Data
% Each model is associated with its data as usual

% Load the data, assigning 'mRNA' and 'RNA_STL1_total_TS3Full' 
% and condition on 'Replica' = 1
STL1_data_rep1 = ...
  STL1_data_rep1.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                        {'mRNA','RNA_STL1_total_TS3Full'},{'Replica',1}); 

% Generate functions for model propensities
STL1_data_rep1 = STL1_data_rep1.formPropensitiesGeneral('STL1_data');

%% Create Second Model and associate to its own data
% Make a copy of the 'STL1_data_rep1' model for Replica 2
STL1_data_rep2 = STL1_data_rep1;

% Load the data, assigning 'mRNA' and 'RNA_STL1_total_TS3Full' 
% and condition on 'Replica' = 2
STL1_data_rep2 = ...
   STL1_data_rep2.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                        {'mRNA','RNA_STL1_total_TS3Full'},{'Replica',2});

%% Set Fitting Options
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display','final','MaxIter',500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(0): Single model
% This is a simple example, where we only fit one model to a single data
% set. First, we create a MultiModel class with just our original model:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

singleModel = SSITMultiModel({STL1_data_rep1},{(1:7)});

% We then copy the original parameters into the MultiModel:
allParsSingle = ([STL1_data_rep1.parameters{:,2}]);

% Next, we run a few rounds of fitting:
for iFit = 1:3
    % Initialize state space:
    singleModel = singleModel.initializeStateSpaces;
    
    % Run seach for MLE:
    allParsSingle = singleModel.maximizeLikelihood(...
        allParsSingle, fitOptions, fitAlgorithm);
    
    % Update Model with new parameters:
    singleModel = singleModel.updateModels(allParsSingle);
end

% We then copy the parameters back into STL1_data_rep1 and STL1_data_rep2 
% so we can reuse them later:
STL1_data_rep1.parameters = singleModel.SSITModels{1}.parameters;
STL1_data_rep2.parameters = singleModel.SSITModels{1}.parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Add a new model+data to an existing multimodel
%   Add a second model + data combination.  In this case, the parameters 
%   of the new model are completely independent of the parameter set for  
%   the first model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combinedModel = singleModel.addModel({STL1_data_rep2},{8:14});
combinedModel = combinedModel.initializeStateSpaces;
allParsCombined = ([STL1_data_rep1.parameters{:,2},...
                   [STL1_data_rep2.parameters{:,2}]]);
allParsCombined = combinedModel.maximizeLikelihood(...
    allParsCombined, fitOptions, fitAlgorithm);
combinedModel = combinedModel.updateModels(allParsCombined);

% Note: This example is shown for illustration purposes only.  Usually, 
% if one is fitting a single model independently, then it is more 
% efficient to fit it separately.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Completely independent parameters
%   Create the combined model in one shot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combinedModelIndependent = ...
    SSITMultiModel({STL1_data_rep1,STL1_data_rep2},{1:7,8:14});
combinedModelIndependent = ...
    combinedModelIndependent.initializeStateSpaces;
allParsIndepdendent = ([STL1_data_rep1.parameters{:,2},...
                       [STL1_data_rep2.parameters{:,2}]]);
allParsIndepdendent = combinedModelIndependent.maximizeLikelihood(...
    allParsIndepdendent, fitOptions, fitAlgorithm);
combinedModelIndependent = ...
    combinedModelIndependent.updateModels(allParsIndepdendent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Completely dependent parameters
%   Use a single set of parameters for both models and data sets. 
%   In the following we make a joint model where both
%   STL1_data_rep1 and STL1_data_rep2 use the parameters [1:7].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combinedModelDependent = SSITMultiModel({STL1_data_rep1,...
                                         STL1_data_rep2},{1:7,1:7});
combinedModelDependent = combinedModelDependent.initializeStateSpaces;
allParsDependent = ([STL1_data_rep1.parameters{:,2}]);
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);
combinedModelDependent = ...
    combinedModelDependent.updateModels(allParsDependent);

% Note: This example is shown for illustration purposes only.  Usually,  
% if one is fitting two replicas of the exact same experiment, then it is
% more efficient to combine the data from both replicas and fit them at 
% the same time, e.g. to combined all replicas into one set, simply load 
% the data as follows:
% STL1_data = ...
%       STL1_data.loadData('data/filtered_data_2M_NaCl_Step.csv',...
%             {'mRNA','RNA_STL1_total_TS3Full'},{'Replica',1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(4): Mixed parameters
%   Sometimes it is desirable to only let some parameters change from
%   condition to condition.  In this example both STL1_data_rep1 and  
%   STL1_data_rep2 use the same parameters [1-4], but parameters [5:7] 
%   are only for STL1_data_rep1 and [8:10] are only for STL1_data_rep2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combinedModelMixed= SSITMultiModel({STL1_data_rep1,STL1_data_rep2},...
                                   {(1:7),[1:4,8:10]});
combinedModelMixed = combinedModelMixed.initializeStateSpaces;
allParsMixed = ([STL1_data_rep1.parameters{:,2},...
                 STL1_data_rep2.parameters{5:7,2}]);
allParsMixed = combinedModelMixed.maximizeLikelihood(allParsMixed, ...
                                               fitOptions, fitAlgorithm);
combinedModelMixed = combinedModelMixed.updateModels(allParsMixed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(5): Constrained parameters
%   It is often helpful to place constraints on parameters, since it can
%   be expected that certain parameters should not change that much from 
%   one experiment to another, while others could be more sensitive to
%   expeimental error.  Here, we will assume that parameters 1-4 are the 
%   same for all cases, and that parameters 5-7 are similar but allowed 
%   to change by small values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constraint = @(x)-sum((x(5:7)-x(8:10)).^2);
combinedModelConstrained = SSITMultiModel({STL1_data_rep1,...
                           STL1_data_rep2},{1:7,[1:4,8:10]},constraint);
combinedModelConstrained = ...
        combinedModelConstrained.initializeStateSpaces;
allParsConstrained = allParsMixed;
allParsConstrained = combinedModelConstrained.maximizeLikelihood(...
                      allParsConstrained, fitOptions, fitAlgorithm);
combinedModelConstrained = ...
        combinedModelConstrained.updateModels(allParsMixed);
