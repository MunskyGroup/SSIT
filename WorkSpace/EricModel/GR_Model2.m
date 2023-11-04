%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations
close all 
clear all
addpath('../CommandLine');

%% Create Model for GR only
ModelGR = SSIT;
ModelGR.species = {'cytGR';'nucGR'};
ModelGR.initialCondition = [20;1];
% ModelGR.propensityFunctions = {'(kcn0+kcn1*IDex)*cytGR';'knc*nucGR'};
% ModelGR.stoichiometry = [-1,1;...
%                          1,-1];
% ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
%     'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01});

% ModelGR.propensityFunctions = {'(kcn0+kcn1*IDex)*cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR'};
% ModelGR.stoichiometry = [-1,1,1,-1,0;...
%                          1,-1,0,0,-1];
ModelGR.propensityFunctions = {'(kcn0+kcn1*IDex)*cytGR';'knc*nucGR';'kg1';'gg1*cytGR'};
ModelGR.stoichiometry = [-1,1,1,-1;...
                         1,-1,0,0];
ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
    'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',0.0000;'gg1',0.00000;'gg2',0.0001});
ModelGR.fspOptions.initApproxSS = true;

ModelGR.dataSet = [];
ModelGR.fittingOptions.modelVarsToFit = [5,6,7,8,9,10,11];
% ModelGR.fittingOptions.modelVarsToFit = [5,6,7,8];

ModelGR100 = ModelGR;
%% Fit using the FSP analysis (100nm Dex)
ModelGR100.inputExpressions = {'IDex','100*exp(-gDex*t)*(t>0)'};
ModelGR100 = ModelGR100.loadData('DUSP1_GR_dataframes/GR_ICC_3hr_Dex_total.csv',...
    {'nucGR','GRNorm'},{'conc','100'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

%%
ModelGR100.solutionScheme = 'FSP';
% ModelGR100.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
%     'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',0.0000;'gg1',0.000000;'gg2',0.0000});
[fspSoln,ModelGR100.fspOptions.bounds] = ModelGR100.solve;
ModelGR100.makePlot(fspSoln,'marginals')
ModelGR100.makeFitPlot

%%
ModelGR100.solutionScheme = 'FSP';
fitOptions = optimset('Display','iter','MaxIter',100);
fitOptions.suppressFSPExpansion = true; 
fitParameters = ModelGR100.maximizeLikelihood([],fitOptions);
ModelGR100.parameters(ModelGR100.fittingOptions.modelVarsToFit,2) = num2cell(fitParameters);
ModelGR100.makeFitPlot

%% Next, build multiple models for different Dex concentrations.
% In these models, the only thing that changes is the Dex concentration,
% which is 100, 10, or 1.
ModelGR10 = ModelGR100;
ModelGR10.inputExpressions = {'IDex','10*exp(-gDex*t)*(t>0)'};
ModelGR10 = ModelGR10.loadData('DUSP1_GR_dataframes/GR_ICC_3hr_Dex_total.csv',...
    {'nucGR','GRNorm'},{'conc','10'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

ModelGR1 = ModelGR100;
ModelGR1.inputExpressions = {'IDex','1*exp(-gDex*t)*(t>0)'};
ModelGR1 = ModelGR1.loadData('DUSP1_GR_dataframes/GR_ICC_3hr_Dex_total.csv',...
    {'nucGR','GRNorm'},{'conc','1'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

GRpars = [ModelGR100.parameters{5:8,2}];

%%  Now to combine all three GR models and fit them.
combinedGRModel = SSITMultiModel({ModelGR1,ModelGR10,ModelGR100},{[1:4],[1:4],[1:4]});
combinedGRModel = combinedGRModel.initializeStateSpaces;
GRpars = combinedGRModel.maximizeLikelihood(...
    GRpars, fitOptions);
combinedGRModel = combinedGRModel.updateModels(GRpars);

%% 
ModelGR1 = combinedGRModel.SSITModels{1};
ModelGR10 = combinedGRModel.SSITModels{10};
ModelGR100 = combinedGRModel.SSITModels{100};

%% Now to extend the fit using an FSP analysis.
ModelGR.solutionScheme = 'FSP';
fitParameters = ModelGR.maximizeLikelihood([],fitOptions);
ModelGR.parameters(ModelGR.fittingOptions.modelVarsToFit,2) = num2cell(fitParameters);
ModelGR.makeFitPlot

%% Load and Associate smFISH Data
% Each model is associated with its data as usual:
Model = SSIT;
Model.species = {'offGene';'onGene';'cytGR';'nucGR';'rna'};
Model.initialCondition = [2;0;24;1;5];
Model.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
    'kcn0+kcn1*IDex)*cytGR';'knc*nucGR';...
    'kr*onGene';'gr*rna'};
Model.stoichiometry = [-1,1,0,0,0,0;...
                         1,-1,0,0,0,0;...
                         0,0,-1,1,0,0;...
                         0,0,1,-1,0,0;...
                         0,0,0,0,1,-1];
Model.inputExpressions = {'IDex','exp(-gDex*t)*(t>0)'};
Model.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
    'kcn0',0.005;'kcn1',0.01;'gDex',0.003;'knc',0.01});
Model.parameters([5,6,7,8],2) = ModelGR.parameters([5,6,7,8],2);

Model.solutionScheme = 'ode';
Model.fspOptions.initApproxSS = true;

Model = Model.loadData('DUSP1_GR_dataframes/DUSP1_3hr_Dex_100nM_total.csv',{'rna','RNA_nuc'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;
Model.solutionScheme = 'ode';
Model.fittingOptions.modelVarsToFit = 1:4;

fitOptions = optimset('Display','iter','MaxIter',100);
fitOptions.suppressFSPExpansion = true; 
fitParameters = Model.maximizeLikelihood([],fitOptions);
Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(fitParameters);
Model.makeFitPlot




%%

%% Create and estimate parameters for combined models
%%      Example 0 -- single model.
% This is a simple example, where we only fit one model to a single data
% set.
singleModel = SSITMultiModel({Model},{[1:7]});
singleModel = singleModel.initializeStateSpaces;
allParsSingle = ([Model.parameters{:,2}]);
allParsSingle = singleModel.maximizeLikelihood(...
    allParsSingle, fitOptions, fitAlgorithm);
singleModel = singleModel.updateModels(allParsSingle);

% Copy the parameters back into Model1 and Model2 so we can reuse them later.
Model.parameters = singleModel.SSITModels{1}.parameters;
ModelGR.parameters = singleModel.SSITModels{1}.parameters;

%%      Example 1 -- Adding new Model+Data to an existing multimodel
% This is how one adds a second model/data combination.  In this case the
% parameters of the new model are completely independent of the parameter
% set for the first model.
combinedModel = singleModel.addModel({ModelGR},{[8:14]});
combinedModel = combinedModel.initializeStateSpaces;
allParsCombined = ([Model.parameters{:,2},[ModelGR.parameters{:,2}]]);
allParsCombined = combinedModel.maximizeLikelihood(...
    allParsCombined, fitOptions, fitAlgorithm);
combinedModel = combinedModel.updateModels(allParsCombined);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting a single model independently, then it is more efficient to fit
% it separately.

%%      Example 2 -- completely independent parameters.
% Here is how we can create the combined model in one shot.
combinedModelIndependent = SSITMultiModel({Model,ModelGR},{[1:7],[8:14]});
combinedModelIndependent = combinedModelIndependent.initializeStateSpaces;
allParsIndepdendent = ([Model.parameters{:,2},[ModelGR.parameters{:,2}]]);
allParsIndepdendent = combinedModelIndependent.maximizeLikelihood(...
    allParsIndepdendent, fitOptions, fitAlgorithm);
combinedModelIndependent = combinedModelIndependent.updateModels(allParsIndepdendent);

%%      Example 3 -- completely dependent parameters.
% Here is an example of how a single set of parameters can be used for both
% models and data sets. In the following we make a joint model where both
% Model1 and Model2 use the parameters [1:7].
combinedModelDependent = SSITMultiModel({Model,ModelGR},{[1:7],[1:7]});
combinedModelDependent = combinedModelDependent.initializeStateSpaces;
allParsDependent = ([Model.parameters{:,2}]);
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);
combinedModelDependent = combinedModelDependent.updateModels(allParsDependent);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting two replicas of the exact same experiment, then it is
% more efficient to combine the data from both replicas and fit them at the
% same time, e.g. to combined all replicas into one set, simply loas the
% data as follows:
% Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});

%%      Example 4 -- mixed parameters.
% Sometimes it is desirable to only let some parameters change from
% condition to condition.  In this example both Model1 and Model2 use the
% same parameters [1-4], but parameters [5:7] are only for Model1 and
% [8:10] are only for Model2.
combinedModelMixed= SSITMultiModel({Model,ModelGR},{[1:7],[1:4,8:10]});
combinedModelMixed = combinedModelMixed.initializeStateSpaces;
allParsMixed = ([Model.parameters{:,2},ModelGR.parameters{5:7,2}]);
allParsMixed = combinedModelMixed.maximizeLikelihood(...
    allParsMixed, fitOptions, fitAlgorithm);

%%      Example 5 -- constrained parameters.
% it is also often helpful to place constraints on parameters since it can
% be expectd that cartain parameters should not change that much from one
% experiment to another.
constraint = @(x)sum((x(5:7)-x(8:10)).^2);
combinedModelConstrained = SSITMultiModel({Model,ModelGR},{[1:7],[1:4,8:10]},constraint);
combinedModelConstrained = combinedModelConstrained.initializeStateSpaces;
allParsConstrained = ([Model.parameters{:,2},ModelGR.parameters{5:7,2}]);
allParsConstrained = combinedModelConstrained.maximizeLikelihood(...
    allParsConstrained, fitOptions, fitAlgorithm);