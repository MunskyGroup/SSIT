%% example_11_ComplexModels_MultiModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Fit multiple models and data sets with shared parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 and 4-state STL1 models from example_1_CreateSSITModels 
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  

% View model summaries:
STL1.summarizeModel
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example script to show how multiple SSIT models and data sets can be fit
%% simultaneously.  This is most useful in situations where:
%   1) The analysis considers different experimental conditions (e.g.,
%      different time points, different inducer concentrations, different
%      genetic mutations).
%   2) Replica-to-replica variations are expected that would result in
%      slightly different parameter combinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make copies of our models:
STL1_multi_1 = STL1;
STL1_4state_multi_1 = STL1_4state;

%% Load and associate smFISH data
%  Each model is associated with its data as usual 
%  (example_8_LoadingandFittingData_DataLoading):

STL1_multi_1 = ...
    STL1_multi_1.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                         {'mRNA','RNA_STL1_total_TS3Full'},...
                         {'Replica',1;'Condition','0.2M_NaCl_Step'});

STL1_4state_multi_1 = ...
   STL1_4state_multi_1.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                               {'mRNA','RNA_STL1_total_TS3Full'},...
                               {'Replica',1;'Condition','0.2M_NaCl_Step'});

% Set up models for FSP solutions:
STL1_multi_1.fspOptions.fspTol = inf;
STL1_4state_multi_1.fspOptions.fspTol = inf;
STL1_multi_1.fittingOptions.modelVarsToFit = 1:7;
STL1_4state_multi_1.fittingOptions.modelVarsToFit = 1:15;
STL1_multi_1.fspOptions.initApproxSS = true;
STL1_4state_multi_1.fspOptions.initApproxSS = true;

% Generate functions for model propensities:
STL1_multi_1 = STL1_multi_1.formPropensitiesGeneral('STL1_multi_FSP');
STL1_4state_multi_1 = ...
    STL1_4state_multi_1.formPropensitiesGeneral('STL1_4state_multi_FSP');

%% Create a second set of models and associate to their own data 
%  In this case, our second set will be associated to Replica 2 data
STL1_multi_2 = STL1_multi_1;
STL1_multi_2 = ...
    STL1_multi_2.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                         {'mRNA','RNA_STL1_total_TS3Full'},...
                         {'Replica',2;'Condition','0.2M_NaCl_Step'});

STL1_4state_multi_2 = STL1_4state_multi_1;
STL1_4state_multi_2 = ...
   STL1_4state_multi_2.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                               {'mRNA','RNA_STL1_total_TS3Full'},...
                               {'Replica',2;'Condition','0.2M_NaCl_Step'});

%% Set Fitting Options
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display','final','MaxIter',500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(0): Solve a single model
%  This is a simple example, where we only fit one model to a single data
%  set, as in example_4_SolveSSITModels_FSP.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STL1_singleModel = SSITMultiModel({STL1_multi_1},{(1:7)});
STL1_4state_singleModel = SSITMultiModel({STL1_4state_multi_1},{(1:15)});

% Copy the original parameters into the MultiModel:
STL1_allParsSingle = ([STL1_multi_1.parameters{:,2}]);
STL1_4state_allParsSingle = ([STL1_4state_multi_1.parameters{:,2}]);

% Run a few rounds of fitting:
for iFit = 1:3
    % Initialize state space:
    STL1_singleModel = STL1_singleModel.initializeStateSpaces;
    STL1_4state_singleModel = ...
        STL1_4state_singleModel.initializeStateSpaces;
    
    % Run seach for MLE:
    STL1_allParsSingle = STL1_singleModel.maximizeLikelihood(...
        STL1_allParsSingle, fitOptions, fitAlgorithm);
    
    STL1_4state_allParsSingle = ...
        STL1_4state_singleModel.maximizeLikelihood(...
        STL1_4state_allParsSingle, fitOptions, fitAlgorithm);
    
    % Update models with new parameters:
    STL1_singleModel = STL1_singleModel.updateModels(STL1_allParsSingle);
    STL1_4state_singleModel = ...
        STL1_4state_singleModel.updateModels(STL1_4state_allParsSingle);
end

% We then copy the parameters back into Model1 and Model2 so we can reuse them
% later:
STL1_multi_1.parameters = STL1_singleModel.SSITModels{1}.parameters;
STL1_multi_2.parameters = STL1_singleModel.SSITModels{1}.parameters;

STL1_4state_multi_1.parameters = ...
    STL1_4state_singleModel.SSITModels{1}.parameters;
STL1_4state_multi_2.parameters = ...
    STL1_4state_singleModel.SSITModels{1}.parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Add a new model+data to existing MultiModel
%  This is how one adds a second model/data combination.  In this case the
%  parameters of the new model are completely independent of the parameter
%  set for the first model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STL1_combinedModel = STL1_singleModel.addModel({STL1_multi_2},{8:14});
STL1_4state_combinedModel = ...
    STL1_4state_singleModel.addModel({STL1_4state_multi_2},{8:14});

combinedModel = combinedModel.initializeStateSpaces;
allParsCombined = ([Model1.parameters{:,2},[Model2.parameters{:,2}]]);
allParsCombined = combinedModel.maximizeLikelihood(...
    allParsCombined, fitOptions, fitAlgorithm);
combinedModel = combinedModel.updateModels(allParsCombined);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting a single model independently, then it is more efficient to fit
% it separately.

%%      Example 2 -- completely independent parameters.
% Here is how we can create the combined model in one shot.
combinedModelIndependent = SSITMultiModel({Model1,Model2},{1:7,8:14});
combinedModelIndependent = combinedModelIndependent.initializeStateSpaces;
allParsIndepdendent = ([Model1.parameters{:,2},[Model2.parameters{:,2}]]);
allParsIndepdendent = combinedModelIndependent.maximizeLikelihood(...
    allParsIndepdendent, fitOptions, fitAlgorithm);
combinedModelIndependent = combinedModelIndependent.updateModels(allParsIndepdendent);

%%      Example 3 -- completely dependent parameters.
% Here is an example of how a single set of parameters can be used for both
% models and data sets. In the following we make a joint model where both
% Model1 and Model2 use the parameters [1:7].
combinedModelDependent = SSITMultiModel({Model1,Model2},{1:7,1:7});
combinedModelDependent = combinedModelDependent.initializeStateSpaces;
allParsDependent = ([Model1.parameters{:,2}]);
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);
combinedModelDependent = combinedModelDependent.updateModels(allParsDependent);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting two replicas of the exact same experiment, then it is
% more efficient to combine the data from both replicas and fit them at the
% same time, e.g. to combined all replicas into one set, simply load the
% data as follows:
% Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'rna','RNA_nuc'});

%%      Example 4 -- mixed parameters.
% Sometimes it is desirable to only let some parameters change from
% condition to condition.  In this example both Model1 and Model2 use the
% same parameters [1-4], but parameters [5:7] are only for Model1 and
% [8:10] are only for Model2.
combinedModelMixed= SSITMultiModel({Model1,Model2},{(1:7),[1:4,8:10]});
combinedModelMixed = combinedModelMixed.initializeStateSpaces;
allParsMixed = ([Model1.parameters{:,2},Model2.parameters{5:7,2}]);
allParsMixed = combinedModelMixed.maximizeLikelihood(...
    allParsMixed, fitOptions, fitAlgorithm);
combinedModelMixed = combinedModelMixed.updateModels(allParsMixed);

%%      Example 5 -- constrained parameters.
% It is often helpful to place constraints on parameters, since it can
% be expected that cartain parameters should not change that much from one
% experiment to another, while others could be more sensitive to
% expeimental error.  Here, we will assume that parameters 1-4 are the same
% for all cases, and that parameters 5-7 are similar but allowed to change
% by small values.
constraint = @(x)-sum((x(5:7)-x(8:10)).^2);
combinedModelConstrained = SSITMultiModel({Model1,Model2},{1:7,[1:4,8:10]},constraint);
combinedModelConstrained = combinedModelConstrained.initializeStateSpaces;
allParsConstrained = allParsMixed;
allParsConstrained = combinedModelConstrained.maximizeLikelihood(...
    allParsConstrained, fitOptions, fitAlgorithm);
combinedModelConstrained = combinedModelConstrained.updateModels(allParsMixed);

%%      Example 6 -- using multimodel to explore batch variations
% In this example, we will use the multimodel to allow parameters to change
% for different replica data sets (e.g., to allow for batch variations, or
% to explore how parameters change under different genetic variations).
% Here, we illustrate a quick means to generate the multimodel starting
% with a single template and a datafile with multiple replicas.

% Sepcify datafile name and species linking rules
DataFileName = '../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv';
LinkedSpecies = {'rna','RNA_nuc'};

% In this case, let's suppose that we only wish to fit the data at times
% before 75 minutes.  We will set the global conditions as:
ConditionsGlobal = {[],[],'TAB.time<=75'};

% We want to split up the replicas to be separate.
ConditionsReplicas = {'TAB.Rep_num==1';...
    'TAB.Rep_num==2'};

% Specify constraints on rep-to-rep parameter variations
Log10Constraints = [0.1,0.1,0.1,0.02,0,0,0];
% Here, we specify that there is an expected 0.1 log10 deviation expected
% in the first three parameters, smaller in the 4th parameter, and no
% deviation at all expected in the last three parameters.

% The full model is created:
CrossValidationModel = SSITMultiModel.createCrossValMultiModel(Model1,DataFileName, ...
                LinkedSpecies,ConditionsGlobal,ConditionsReplicas, ...
                Log10Constraints);
CrossValidationModel = CrossValidationModel.initializeStateSpaces;

% Now to run the model fitting routines
crossValPars = CrossValidationModel.parameters;
crossValPars = CrossValidationModel.maximizeLikelihood(...
    crossValPars, fitOptions, fitAlgorithm);
CrossValidationModel = CrossValidationModel.updateModels(crossValPars);
CrossValidationModel.parameters = crossValPars;

% Make a figure to explore how much the parameters changed between replicas.
fignum = 12; useRelative = true;
CrossValidationModel.compareParameters(fignum,useRelative);