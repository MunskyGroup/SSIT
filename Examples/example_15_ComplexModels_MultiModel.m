%% SSIT/Examples/example_15_ComplexModels_MultiModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Complex models
%   * Fit multiple models and data sets with shared parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the solved 4-state STL1 model from example_4_SolveSSITModels_FSP 
%clear
%close all

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

%% Set Fitting Options:
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display','final','MaxIter',200); % small # for demo

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

% Select which models to include in SSITMultiModel:
Models_ind = {STL1_4state_multi_1_ind, STL1_4state_multi_2_ind};

%% Define how parameters are assigned to sub-models by their indices.  
% In this example, the parameters are completely independent:
ParsIndices_ind = {[1:13], [14:26]};

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

% Select which models to include in SSITMultiModel:
Models_dep = {STL1_4state_multi_1_dep, STL1_4state_multi_2_dep};

%% Define how parameters are assigned to sub-models by their indices.  
% In this example, the parameters are completely dependent:
ParsIndices_dep = {[1:13], [1:13]};

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
% STL1_4state_multi_2 share parameters [1:11], but each model also has its
% own [12:13]. (This was shown in 'example_11_...CrossValidation.m')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make copies of our multi models:
STL1_4state_multi_1_mix = STL1_4state_multi_1;
STL1_4state_multi_2_mix = STL1_4state_multi_2;

% Select which models to include in SSITMultiModel:
Models_mix = {STL1_4state_multi_1_mix, STL1_4state_multi_2_mix};

%% Define how parameters are assigned to sub-models by their indices  
% In this example, the first 11 parameters are shared, and each model has
% two of its own parameters which must be stored in separate indices:
ParsIndices_mix = {[1:11,12:13],[1:11,14:15]};

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
%% Ex(5): Different models, same data
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
% Model 1 has one shared parameter plus 6 of its own (7 total):
STL1_multi.fittingOptions.modelVarsToFit = [4,[1:3,5:7]];

% Model 2 has one shared parameter plus 12 of its own (13 total):
STL1_4state_multi.fittingOptions.modelVarsToFit = [9,[1:8,10:13]];

% Select which models to include in SSITMultiModel:
Models_diff = {STL1_multi, STL1_4state_multi};

%% Define how parameters are assigned to sub-models by their indices  
ParsIndices_diff = {[1,2:7], [1,8:19]};

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