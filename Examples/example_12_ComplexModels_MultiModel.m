%% example_12_ComplexModels_MultiModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Fit multiple models and data sets with shared parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 and 4-state STL1 models from example_1_CreateSSITModels 
%clear
%close all
addpath(genpath('../'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE
% example_10_LoadingandFittingData_MHA

%% Load pre-computed model fit
% load('example_9_LoadingandFittingData_MLE.mat')

% View model summariy:
STL1_4state_MLE.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example script to show how multiple SSIT models and data sets can be fit
%% simultaneously.  This is most useful in situations where:
%   1) The analysis considers different experimental conditions (e.g.,
%      different time points, different inducer concentrations, different
%      genetic mutations).
%   2) Replica-to-replica variations are expected that would result in
%      slightly different parameter combinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make a copy of our model:
STL1_4state_multi_1 = STL1_4state_MLE;

%% Load and associate smFISH data
%  Associate the data with an SSIT model data as usual 
%  (example_8_LoadingandFittingData_DataLoading):

STL1_4state_multi_1 = ...
   STL1_4state_multi_1.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                               {'mRNA','RNA_STL1_total_TS3Full'},...
                               {'Replica',1;'Condition','0.2M_NaCl_Step'});


%% Create a second model and associate it to its own data 
%  In this case, the second set will be associated to Replica 2 data

STL1_4state_multi_2 = STL1_4state_multi_1;
STL1_4state_multi_2 = ...
   STL1_4state_multi_2.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                               {'mRNA','RNA_STL1_total_TS3Full'},...
                               {'Replica',2;'Condition','0.2M_NaCl_Step'});

%% Set Fitting Options
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display','final','MaxIter',200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(0): Single model
% This is a simple example, where we only fit one model to a single data
% set. First, we create a MultiModel class with just our original model:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npars = size(STL1_4state_multi_1.parameters,1);
STL1_4state_multi_1.fittingOptions.modelVarsToFit = 1:npars;   
% or your intended subset, all ≤ npar
singleModel = SSITMultiModel({STL1_4state_multi_1},{1:npars});

% We then copy the original parameters into the MultiModel:
allParsSingle = ([STL1_4state_multi_1.parameters{:,2}]);

% Next, we run a few rounds of fitting:
for iFit = 1:2
    % Initialize state space:
    singleModel = singleModel.initializeStateSpaces;
    
    % Run seach for MLE:
    allParsSingle = singleModel.maximizeLikelihood(...
        allParsSingle, fitOptions, fitAlgorithm);
    
    % Update Model with new parameters:
    singleModel = singleModel.updateModels(allParsSingle);
end

% We then copy the parameters back into Model1 and Model2 so we can reuse 
% them later:
STL1_4state_multi_1.parameters = singleModel.SSITModels{1}.parameters;
STL1_4state_multi_2.parameters = singleModel.SSITModels{1}.parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Adding new model+data to an existing MultiModel
% This is how one adds a second model/data combination.  In this case the
% parameters of the new model are completely independent of the parameter
% set for the first model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npars2 = size(STL1_4state_multi_2.parameters,1);
STL1_4state_multi_2.fittingOptions.modelVarsToFit = 1:npars2;   
% or your intended subset, all ≤ npar

combinedModel = singleModel.addModel({STL1_4state_multi_2},{1:npars2});
combinedModel = combinedModel.initializeStateSpaces;
allParsCombined = ([STL1_4state_multi_1.parameters{:,2},...
                   [STL1_4state_multi_2.parameters{:,2}]]);
allParsCombined = combinedModel.maximizeLikelihood(...
    allParsCombined, fitOptions, fitAlgorithm);
combinedModel = combinedModel.updateModels(allParsCombined);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting a single model independently, then it is more efficient to 
% fit it separately.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Completely independent parameters
%  Here is how we can create the combined model in one shot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combinedModelIndependent = SSITMultiModel({STL1_4state_multi_1,...
                                        STL1_4state_multi_2},{1:15,16:30});
combinedModelIndependent = combinedModelIndependent.initializeStateSpaces;
allParsIndepdendent = ([STL1_4state_multi_1.parameters{:,2},...
                       [STL1_4state_multi_2.parameters{:,2}]]);
allParsIndepdendent = combinedModelIndependent.maximizeLikelihood(...
    allParsIndepdendent, fitOptions, fitAlgorithm);
combinedModelIndependent = ...
    combinedModelIndependent.updateModels(allParsIndepdendent);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Completely dependent parameters
%  Here is an example of how a single set of parameters can be used for 
%  both models and data sets. In the following we make a joint model where 
%  both STL1_4state_multi_1 and STL1_4state_multi_2 use the parameters 
% [1:15].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combinedModelDependent = SSITMultiModel({STL1_4state_multi_1,...
                                         STL1_4state_multi_2},{1:15,1:15});
combinedModelDependent = combinedModelDependent.initializeStateSpaces;
allParsDependent = ([STL1_4state_multi_1.parameters{:,2}]);
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);
combinedModelDependent = ...
    combinedModelDependent.updateModels(allParsDependent);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting two replicas of the exact same experiment, then it is
% more efficient to combine the data from both replicas and fit them at the
% same time, e.g. to combined all replicas into one set, simply load the
% data without the replica filtering condition:
% STL1_4state_data = ...
%     STL1_4state_data.loadData('data/filtered_data_2M_NaCl_Step.csv',...
%                              {'mRNA','RNA_STL1_total_TS3Full'},...
%                              {'Condition','0.2M_NaCl_Step'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(4): Mixed parameters
% Sometimes it is desirable to only let some parameters change from
% condition to condition.  In this example both STL1_4state_multi_1 and  
% STL1_4state_multi_2 use the same parameters [1-11], but parameters 
% [12:13] are only for STL1_4state_multi_1 and [14:15] are only for 
% STL1_4state_multi_2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model 1 contributes shared [1:11] plus its own [12:13]
STL1_4state_multi_1.fittingOptions.modelVarsToFit = [1:11, 12:13];

% Model 2 contributes shared [1:11] plus its own [14:15]
STL1_4state_multi_2.fittingOptions.modelVarsToFit = [1:11, 14:15];

combinedModelMixed = SSITMultiModel( ...
    {STL1_4state_multi_1, STL1_4state_multi_2}, ...
    { [1:11, 12:13],       [1:11, 14:15] } );
combinedModelMixed = combinedModelMixed.initializeStateSpaces;
allParsMixed = ([STL1_4state_multi_1.parameters{:,2},...
                 STL1_4state_multi_2.parameters{7:15,2}]);
allParsMixed = combinedModelMixed.maximizeLikelihood(...
    allParsMixed, fitOptions, fitAlgorithm);
combinedModelMixed = combinedModelMixed.updateModels(allParsMixed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(5): Constrained parameters
% It is often helpful to place constraints on parameters, since it can
% be expected that certain parameters should not change that much from one
% experiment to another, while others could be more sensitive to
% expeimental error.  Here, we will assume that parameters 1-6 are the same
% for all cases, and that parameters 7-15 are similar but allowed to change
% by small values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constraint = @(x)-sum((x(7:15)-x(8:16)).^2);
combinedModelConstrained = SSITMultiModel({STL1_4state_multi_1,...
    STL1_4state_multi_2},{1:16,[1:4,8:16]},constraint);
combinedModelConstrained = combinedModelConstrained.initializeStateSpaces;
allParsConstrained = allParsMixed;
allParsConstrained = combinedModelConstrained.maximizeLikelihood(...
    allParsConstrained, fitOptions, fitAlgorithm);
combinedModelConstrained = ...
    combinedModelConstrained.updateModels(allParsMixed);