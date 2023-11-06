%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
clear all
clc
%% Define SSIT Model
  Model1 = SSIT;
  Model1.species = {'x1';'x2'};
  Model1.initialCondition = [0;0];
  Model1.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
  Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
  Model1.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
  Model1.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01;...
                       'a1',0.4;'r1',0.04;'r2',0.1});
  Model1.fspOptions.initApproxSS = true;

%% Solve the model using the FSP
Model1.solutionScheme = 'FSP';
Model1.fspOptions.fspTol = 1e-4;
% Model.fspOptions.bounds(4) = 300;
[fspSoln,Model1.fspOptions.bounds] = Model1.solve;

%% Load and Fit smFISH Data
% Model = Model.loadData('../../DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc';'x3','GR_nuc'},...
%     {'Rep_num','1'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;
Model1 = Model1.loadData('../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'},...
    {'Rep_num','1'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

Model1.fspOptions.fspTol = inf;
Model1.fittingOptions.modelVarsToFit = 1:7;
fitOptions = optimset('Display','iter','MaxIter',50);
%% fit the first model to first data set
Model1.parameters(Model1.fittingOptions.modelVarsToFit,2) = num2cell(Model1.maximizeLikelihood([],fitOptions));
Model1.makeFitPlot;

%% Create Second Model
Model2 = Model1;

% And change accordingly.
Model2 = Model2.loadData('../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'},...
    {'Rep_num','2'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;

%% Create a vector of shared parameters
allPars = log10([Model1.parameters{:,2},Model2.parameters{:,2}]);

%% Create a shared objective function to minimize
objTotal = @(x)-Model1.computeLikelihood(10.^x(1:7),fspSoln.stateSpace) + ...
    -Model2.computeLikelihood(10.^x(8:14),fspSoln.stateSpace);
objTotal(allPars)

%% Search parameter space to fit both models to both data sets at the same time.
allPars = fminsearch(objTotal,allPars,fitOptions);

%% Update and make plots of model results
Model1.parameters(:,2)  = num2cell(10.^allPars(1:7));
Model2.parameters(:,2)  = num2cell(10.^allPars(8:14));
Model1.makeFitPlot
Model2.makeFitPlot

%% Add constraints on parameters
% option 1 -- both models use exact same parameters
% allParsOpt1 = allPars(1:7);
objTotal = @(x)-Model1.computeLikelihood(10.^x(1:7),fspSoln.stateSpace) + ...
    -Model2.computeLikelihood(10.^x(1:7),fspSoln.stateSpace);
objTotal(allParsOpt1)
allParsOpt1 = fminsearch(objTotal,allParsOpt1,fitOptions);
Model1.parameters(:,2)  = num2cell(10.^allParsOpt1(1:7));
Model2.parameters(:,2)  = num2cell(10.^allParsOpt1(1:7));
Model1.makeFitPlot
Model2.makeFitPlot

%% Add constraints on parameters
% option 2 -- parameters can change, but only if they need to
% allParsOpt2 = allPars(1:14);
allParsOpt2 = [allParsOpt1,allParsOpt1];
lambda = 100;
objTotal2 = @(x)-Model1.computeLikelihood(10.^x(1:7),fspSoln.stateSpace) + ...
    -Model2.computeLikelihood(10.^x(8:14),fspSoln.stateSpace)+...
    lambda*sum(abs(x(1:7)-x(8:14)));
objTotal2(allParsOpt2)
allParsOpt2 = fminsearch(objTotal2,allParsOpt2,fitOptions);
Model1.parameters(:,2)  = num2cell(10.^allParsOpt2(1:7));
Model2.parameters(:,2)  = num2cell(10.^allParsOpt2(8:14));
Model1.makeFitPlot
Model2.makeFitPlot
[10.^allParsOpt2(1:7)',10.^allParsOpt2(8:14)']

%% Add constraints on parameters
% option 3 -- only some parameters can change, and only if they need to
% allParsOpt2 = allPars(1:14);
allParsOpt3 = [allParsOpt1(1:7),allParsOpt1(5:7)];
lambda = 100;
objTotal3 = @(x)-Model1.computeLikelihood(10.^x(1:7),fspSoln.stateSpace) + ...
    -Model2.computeLikelihood(10.^x([1:4,8:10]),fspSoln.stateSpace)+...
    lambda*sum(abs(x(5:7)-x(8:10)));
objTotal3(allParsOpt3)
allParsOpt3 = fminsearch(objTotal3,allParsOpt3,fitOptions);
Model1.parameters(:,2)  = num2cell(10.^allParsOpt3(1:7));
Model2.parameters(:,2)  = num2cell(10.^allParsOpt3([1:4,8:10]));
Model1.makeFitPlot
Model2.makeFitPlot
[10.^allParsOpt3(1:7)',10.^allParsOpt3([1:4,8:10])']