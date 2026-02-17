%% SSIT/Examples/example_15_ComplexModels_MultiModel

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

%% Load pre-solved model 
% load('example_4_SolveSSITModels_FSP.mat')

% View model summariy:
STL1_4state_FSP.summarizeModel

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
STL1_4state_multi_1 = STL1_4state_FSP;

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

%% Set Fitting Options:
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display','final','MaxIter',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(0): Single model
% This is a simple example, where we only fit one model to a single data
% set. First, we create a MultiModel class with just our original model:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make copies of our multi models:
STL1_4state_multi_1_single = STL1_4state_multi_1;
STL1_4state_multi_2_single = STL1_4state_multi_2;

npars = size(STL1_4state_multi_1_single.parameters,1);
STL1_4state_multi_1_single.fittingOptions.modelVarsToFit = 1:npars;   
% or your intended subset, all ≤ npar
singleModel = SSITMultiModel({STL1_4state_multi_1_single},{1:npars});

% We then copy the original parameters into the MultiModel:
allParsSingle = ([STL1_4state_multi_1_single.parameters{:,2}]);

% Next, we run a few rounds of fitting:
for iFit = 1:2
    % Initialize state space:
    singleModel = singleModel.initializeStateSpaces;
    
    % Run seach for MLE:
    allParsSingle = singleModel.maximizeLikelihood(allParsSingle,...
                                                 fitOptions, fitAlgorithm);
    
    % Update Model with new parameters:
    singleModel = singleModel.updateModels(allParsSingle);
end

% We then copy the parameters back into Model1 and Model2 so we can reuse 
% them later:
STL1_4state_multi_1_single.parameters = singleModel.SSITModels{1}.parameters;
STL1_4state_multi_2_single.parameters = singleModel.SSITModels{1}.parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Adding new model+data to an existing MultiModel
% This is how one adds a second model/data combination.  In this case the
% parameters of the new model are completely independent of the parameter
% set for the first model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make copies of our multi models:
STL1_4state_multi_1_new = STL1_4state_multi_1;
STL1_4state_multi_2_new = STL1_4state_multi_2;

npars2 = size(STL1_4state_multi_2_new.parameters,1);
STL1_4state_multi_2_new.fittingOptions.modelVarsToFit = 1:npars2;   
% or your intended subset, all ≤ npar

combinedModel = singleModel.addModel({STL1_4state_multi_2_new},{1:npars2});
combinedModel = combinedModel.initializeStateSpaces;
allParsCombined = ([STL1_4state_multi_1_new.parameters{:,2},...
                   [STL1_4state_multi_2_new.parameters{:,2}]]);
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

% Make copies of our multi models:
STL1_4state_multi_1_ind = STL1_4state_multi_1;
STL1_4state_multi_2_ind = STL1_4state_multi_2;

%% Specify how many model parameters will be fit
% Model 1 has 0 shared parameters plus 18 of its own (18 total):
STL1_4state_multi_1_ind.fittingOptions.modelVarsToFit = [1:18];

% Model 2 has 0 shared parameters plus 18 of its own (18 total):
STL1_4state_multi_2_ind.fittingOptions.modelVarsToFit = [1:18];

% Select which models to include in SSITMultiModel:
Models_ind = {STL1_4state_multi_1_ind, STL1_4state_multi_2_ind};

% Define how parameters are assigned to sub-models by their indices.  
% In this example, the parameters are completely independent:
ParsIndices_ind = {[1:18], [19:36]};

% Combine models into one "MultiModel", specify parameters, and initialize:
combinedModelIndependent = SSITMultiModel(Models_ind,ParsIndices_ind);
combinedModelIndependent = combinedModelIndependent.initializeStateSpaces;

% Store parameters for later updating:
allParsIndependent = ([STL1_4state_multi_1_ind.parameters{:,2},...
                       [STL1_4state_multi_2_ind.parameters{:,2}]]);

% Fit parameters using maximum likelihood estimation:
allParsIndependent = combinedModelIndependent.maximizeLikelihood(...
    allParsIndependent, fitOptions, fitAlgorithm);

% Update model parameters and plot results:
combinedModelIndependent = ...
    combinedModelIndependent.updateModels(allParsIndependent,true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Completely dependent parameters
%  Here is an example of how a single set of parameters can be used for 
%  both models and data sets. In the following we make a joint model where 
%  both STL1_4state_multi_1 and STL1_4state_multi_2 use the parameters 
% [1:15].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make copies of our multi models:
STL1_4state_multi_1_dep = STL1_4state_multi_1;
STL1_4state_multi_2_dep = STL1_4state_multi_2;

%% Specify how many model parameters will be fit
% Model 1 has 18 shared parameters plus 0 of its own (18 total): 
STL1_4state_multi_1_dep.fittingOptions.modelVarsToFit = [1:18];

% Model 2 has 18 shared parameters plus 0 of its own (18 total):
STL1_4state_multi_2_dep.fittingOptions.modelVarsToFit = [1:18];

% Select which models to include in SSITMultiModel:
Models_dep = {STL1_4state_multi_1_dep, STL1_4state_multi_2_dep};

% Define how parameters are assigned to sub-models by their indices.  
% In this example, the parameters are completely dependent:
ParsIndices_dep = {[1:18], [1:18]};

% Combine models into one "MultiModel", specify parameters, and initialize:
combinedModelDependent = SSITMultiModel(Models_dep, ParsIndices_dep);
combinedModelDependent = combinedModelDependent.initializeStateSpaces;

% Store parameters for later updating:
allParsDependent = ([STL1_4state_multi_1_dep.parameters{:,2}]);

% Fit parameters using maximum likelihood estimation:
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);

% Update model parameters and plot results:
combinedModelDependent = ...
    combinedModelDependent.updateModels(allParsDependent,true);

%% Note: This example is shown for illustration purposes only.  
% Usually, if one is fitting two replicas of the exact same experiment, 
% then it is more efficient to combine the data from both replicas and fit 
% them at the same time, e.g. to combined all replicas into one set, simply 
% load the data without the replica filtering condition:
% STL1_4state_data = ...
%     STL1_4state_data.loadData('data/filtered_data_2M_NaCl_Step.csv',...
%                              {'mRNA','RNA_STL1_total_TS3Full'},...
%                              {'Condition','0.2M_NaCl_Step'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(4): Mixed parameters
% Sometimes it is desirable to only let some parameters change from
% condition to condition.  In this example, both STL1_4state_multi_1 and  
% STL1_4state_multi_2 share parameters [14:18], but each model also has
% its own [1:13].  (This was shown in 'example_11_...CrossValidation.m')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make copies of our multi models:
STL1_4state_multi_1_mix = STL1_4state_multi_1;
STL1_4state_multi_2_mix = STL1_4state_multi_2;

%% Specify how many model parameters will be fit
% Model 1 has 5 shared parameters plus 13 of its own (18 total):
STL1_4state_multi_1_mix.fittingOptions.modelVarsToFit = [14:18,1:13];

% Model 2 has 5 shared parameters plus 13 of its own (18 total):
STL1_4state_multi_2_mix.fittingOptions.modelVarsToFit = [14:18,1:13];

% Select which models to include in SSITMultiModel:
Models_mix = {STL1_4state_multi_1_mix, STL1_4state_multi_2_mix};

%% Define how parameters are assigned to sub-models by their indices  
% In this example, the first 13 parameters are shared, and each model has 5
% of its own parameters which must be stored in separate indices:
ParsIndices_mix = {[1:13,14:18], [1:13,19:23]};

% Combine models into one "MultiModel", specify parameters, and initialize:
combinedModelMixed = SSITMultiModel(Models_mix, ParsIndices_mix);
combinedModelMixed = combinedModelMixed.initializeStateSpaces;

% Store parameters for later updating:
allParsMixed = ([STL1_4state_multi_1_mix.parameters{:,2},...
                 STL1_4state_multi_2_mix.parameters{:,2}]);

% Fit parameters using maximum likelihood estimation:
allParsMixed = combinedModelMixed.maximizeLikelihood(...
                allParsMixed, fitOptions, fitAlgorithm);

% Update model parameters and plot results:
combinedModelMixed = combinedModelMixed.updateModels(allParsMixed, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(5): Constrained parameters
% It is often helpful to place constraints on parameters, since it can
% be expected that certain parameters should not change that much from one
% experiment to another, while others could be more sensitive to
% expeimental error.  Here, we will assume that parameters 14-18 are the 
% same for both cases, and that parameters 1-13 are similar but allowed to 
% change by small values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make copies of our multi models:
STL1_4state_multi_1_con = STL1_4state_multi_1;
STL1_4state_multi_2_con = STL1_4state_multi_2;

% Each STL1_4state_multi_* model has 18 fit parameters locally
STL1_4state_multi_1_con.fittingOptions.modelVarsToFit = 1:18;
STL1_4state_multi_2_con.fittingOptions.modelVarsToFit = 1:18;

%% Specify the constraint:
sigma = 1.0;
constraint = @(x) -(1/(2*sigma^2)) * sum( (x(1:13) - x(19:31)).^2 );

%% Define how parameters are assigned to sub-models by their indices  
parIdx = { 1:18, [14:18, 19:31] };

% Combine models into one "MultiModel", specify parameters, and initialize:
combinedModelConstrained = SSITMultiModel({STL1_4state_multi_1_con,...
                           STL1_4state_multi_2_con}, parIdx, constraint);

combinedModelConstrained = combinedModelConstrained.initializeStateSpaces();

x0 = allParsMixed;      % 1..15 exist
x0(16:24) = x0(7:15);   % seed the second block

% Fit parameters using maximum likelihood estimation:
xOpt = combinedModelConstrained.maximizeLikelihood(x0, fitOptions,...
                                                   fitAlgorithm);

% Update model parameters and plot results:
combinedModelConstrained = combinedModelConstrained.updateModels(xOpt,true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(6): Different models, same data
% In this example, both the simplified STL1 model "STL1_multi" and the
% "STL1_4state_multi" model share the parameter "dr", [4] and [9], 
% respectively, but parameters [1:3,5:7] are only for Model 1 (STL1) and 
% [1:8,10:18] are only for Model 2 (STL1_4state).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make copies of our models:
STL1_multi = STL1_FSP;
STL1_4state_multi = STL1_4state_multi_1;

% Load and associate data:
STL1_multi = STL1_multi.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                               {'mRNA','RNA_STL1_total_TS3Full'},...
                               {'Replica',1;'Condition','0.2M_NaCl_Step'});
STL1_4state_multi = ...
   STL1_4state_multi.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                             {'mRNA','RNA_STL1_total_TS3Full'},...
                             {'Replica',1;'Condition','0.2M_NaCl_Step'});

%% Specify how many model parameters will be fit
% Model 1 has 1 shared parameters plus 6 of its own (7 total):
STL1_multi.fittingOptions.modelVarsToFit = [4,[1:3,5:7]];

% Model 2 has 1 shared parameters plus 14 of its own (15 total):
STL1_4state_multi.fittingOptions.modelVarsToFit = [9,[1:8,10:18]];

% Select which models to include in SSITMultiModel:
Models_diff = {STL1_multi, STL1_4state_multi};

%% Define how parameters are assigned to sub-models by their indices  
ParsIndices_diff = {[1,2:7], [1,8:24]};

% Combine models into one "MultiModel", specify parameters, and initialize:
combinedModeldiff = SSITMultiModel(Models_diff, ParsIndices_diff);
combinedModeldiff = combinedModeldiff.initializeStateSpaces;

% Store parameters for later updating:
allParsDiff = ([STL1_multi.parameters{:,2},...
                STL1_4state_multi.parameters{:,2}]);

% Fit parameters using maximum likelihood estimation:
allParsDiff = combinedModeldiff.maximizeLikelihood(allParsDiff,...
                                                 fitOptions, fitAlgorithm);

% Update model parameters and plot results:
combinedModeldiff = combinedModeldiff.updateModels(allParsDiff, true);