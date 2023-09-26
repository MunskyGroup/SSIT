%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations

clear all
clc

%% Define SSIT Model
% SSIT models are defined as usual:
Model1 = SSIT;
Model1.species = {'x1';'x2'};
Model1.initialCondition = [0;0];
Model1.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model1.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
Model1.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01;...
    'a1',0.4;'r1',0.04;'r2',0.1});
Model1.fspOptions.initApproxSS = true;

%% Load and Associate smFISH Data
% Each model is associated with its data as usual:
Model1 = Model1.loadData('../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'},...
    {'Rep_num','1'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

Model1.fspOptions.fspTol = inf;
Model1.fittingOptions.modelVarsToFit = 1:7;

%% Create Second Model and associate to its own data
Model2 = Model1;
Model2 = Model2.loadData('../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'},...
    {'Rep_num','2'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting two replicas of the same experiment, then it is
% more efficient to combine the data from both replicas and fit them at the
% same time, e.g. to combined all replicas into one set, simply loas the
% data as follows:
% Model = Model.loadData('../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});

%% Fitting options
fitOptions = optimset('Display','iter','MaxIter',100);
fitAlgorithm = 'fminsearch';

%% Create and estimate parameters for combined models
%%      Example 0 -- single model.
singleModel = SSITMultiModel({Model1},{[1:7]});
singleModel = singleModel.initializeStateSpaces;
allParsSingle = ([Model1.parameters{:,2}]);
allParsSingle = singleModel.maximizeLikelihood(...
    allParsSingle, fitOptions, fitAlgorithm);
singleModel.updateModels(allParsSingle);

%%      Example 1 -- Adding new Model+Data to an existing multimodel
% This is how one adds a second model/data combination.  In this case the
% parameters of the new model are completely independent of the parameteer
% set for the first model.
combinedModel = singleModel.addModel({Model2},{[8:14]});
combinedModel = combinedModel.initializeStateSpaces;
allParsCombined = ([Model1.parameters{:,2},[Model2.parameters{:,2}]]);
allParsCombined = combinedModel.maximizeLikelihood(...
    allParsCombined, fitOptions, fitAlgorithm);
combinedModel.updateModels(allParsCombined);

%%      Example 2 -- completely independent parameters.
% Here is how we can create the combined model in one shot.
combinedModelIndependent = SSITMultiModel({Model1,Model2},{[1:7],[8:14]});
combinedModelIndependent = combinedModelIndependent.initializeStateSpaces;
allParsIndepdendent = ([Model1.parameters{:,2},[Model2.parameters{:,2}]]);
allParsIndepdendent = combinedModelIndependent.maximizeLikelihood(...
    allParsIndepdendent, fitOptions, fitAlgorithm);
combinedModelIndependent.updateModels(allParsIndepdendent);

%%      Example 3 -- completely dependent parameters.
% Here is an example of how a single set of parameters can be used for both
% models and data sets.
combinedModelDependent = SSITMultiModel({Model1,Model2},{[1:7],[1:7]});
combinedModelDependent = combinedModelDependent.initializeStateSpaces;
allParsDependent = ([Model1.parameters{:,2}]);
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);
combinedModelDependent.updateModels(allParsDependent);

%%      Example 4 -- mixed parameters.
% Sometimes it is desirable to only let some parameters change from
% condition to condition.
combinedModelMixed= SSITMultiModel({Model1,Model2},{[1:7],[1:4,8:10]});
combinedModelMixed = combinedModelMixed.initializeStateSpaces;
allParsMixed = ([Model1.parameters{:,2},Model2.parameters{5:7,2}]);
allParsMixed = combinedModelMixed.maximizeLikelihood(...
    allParsMixed, fitOptions, fitAlgorithm);

%%      Example 5 -- constrained parameters.
% it is also often helpful to place constraints on parameters since it can
% be expectd that cartain parameters should not change that much from one
% experiment to another.
constraint = @(x)sum((x(5:7)-x(8:10)).^2);
combinedModelConstrained = SSITMultiModel({Model1,Model2},{[1:7],[1:4,8:10]},constraint);
combinedModelConstrained = combinedModelConstrained.initializeStateSpaces;
allParsConstrained = ([Model1.parameters{:,2},Model2.parameters{5:7,2}]);
allParsConstrained = combinedModelConstrained.maximizeLikelihood(...
    allParsConstrained, fitOptions, fitAlgorithm);