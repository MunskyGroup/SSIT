%% example_12_ComplexModels_MultiModel

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
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE
% example_10_LoadingandFittingData_MHA

% View model summariy:
STL1_4state_MH_it.summarizeModel

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
STL1_4state_multi_1 = STL1_4state_MH_it;

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
fitOptions = optimset('Display','final','MaxIter',500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Completely dependent parameters
%  Here is an example of how a single set of parameters can be used for 
%  both models and data sets. In the following we make a joint model where 
%  both Model1 and Model2 use the parameters [1:15].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combinedModel_Dependent = ...
    SSITMultiModel({STL1_4state_multi_1,STL1_4state_multi_2},{1:15,1:15});
combinedModel_Dependent = combinedModel_Dependent.initializeStateSpaces;
allPars_Dependent = ([STL1_4state_multi_1.parameters{:,2}]);
allPars_Dependent = ...
    combinedModel_Dependent.maximizeLikelihood(allPars_Dependent,...
                                               fitOptions,fitAlgorithm);
combinedModel_Dependent = ...
    combinedModel_Dependent.updateModels(allPars_Dependent);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting two replicas of the exact same experiment, then it is
% more efficient to combine the data from both replicas and fit them at the
% same time, e.g. to combined all replicas into one set, simply load the
% data without setting a condition to filter by replica, e.g. -

% STL1_4state_multi_1 = ...
%   STL1_4state_multi_1.loadData('data/filtered_data_2M_NaCl_Step.csv',...
%                               {'mRNA','RNA_STL1_total_TS3Full'},...
%                               {'Condition','0.2M_NaCl_Step'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Completely independent parameters
%  This is how one adds a second model/data combination.  In this case the
%  parameters of the new model are completely independent of the parameter
%  set for the first model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set propensity functions:
% STL1_4state_multi_2.propensityFunctions = ...
%     {'kr1 * s1';'k12 * s1';'k21 * s2 * IHog';...
%      'kr2 * s2';'k23 * s2';'k32 * s3';...
%      'kr3 * s3';'k34 * s3';'k43 * s4';...
%      'kr4 * s4';'dr * mRNA'}; 
% 
% % Create separate parameters for Model2:
STL1_4state_multi_2.parameters = ...
    ({'j12',0; 'j23',0; 'j34',0;...
      'j21',0; 'j32',0; 'j43',0; ...
      'jr1',0; 'jr2',0; 'jr3',0; 'jr4',0; ...
      'gr',0; 'b0',0; 'b1',0; 'm1',0; 'm2',0;...
      'k12',0.2; 'k23',0.2; 'k34',0.5;...
      'k21',1e-03; 'k32',5e-04; 'k43',5e-04; ...
      'kr1',0.1; 'kr2',0.1; 'kr3',20; 'kr4',0.1; ...
      'dr',0.3; 'a0',1e-03; 'a1',1e-03; 'r1',1e-03; 'r2',1e-03});
% 
% STL1_4state_multi_2.parameters{16:30} = ...
%     STL1_4state_multi_1.parameters{1:15};

% F = F.addParameter({'kr',0.1})

% Here is how we can create the combined model in one shot.
combinedModel_Independent = ...
    SSITMultiModel({STL1_4state_multi_1,STL1_4state_multi_2},{1:15,16:30});
combinedModel_Independent = ...
    combinedModel_Independent.initializeStateSpaces;
allPars_Independent = ([STL1_4state_multi_1.parameters{:,2},...
                        [STL1_4state_multi_2.parameters{:,2}]]);
allPars_Independent = combinedModel_Independent.maximizeLikelihood(...
    allPars_Independent, fitOptions, fitAlgorithm);
combinedModel_Independent = ...
    combinedModel_Independent.updateModels(allPars_Independent);

%%      Example 5 -- constrained parameters.
% It is often helpful to place constraints on parameters, since it can
% be expected that certain parameters should not change that much from one
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