%% example_6_LoadingandFittingData_MLE
% Example script to demonstrate how to load and fit
% experimental data using maximum likelihood estimates (MLEs)
clear
close all
addpath(genpath('../src'));

%% Preliminaries
% Load our STL1 Model described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
example_2_SolveSSITModels_FSP
Model.summarizeModel
STL1Model.summarizeModel

% Set fitOptions, with the maximum allowable number of iterations to fit
fitOptions = optimset('Display','iter','MaxIter',500);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
Modelpars = cell2mat(Model.parameters(1:4,2));
STL1pars = cell2mat(STL1Model.parameters(1:8,2));

%% Load experimental data, matching the species name of the model ('mRNA') 
% to the appropriate column in the data file ('rna')
STL1Real = STL1Model.loadData('data/STL1.csv',{'mRNA','rna'});

% These plots are unnecessary, as the model parameters have not been fit to
% the data yet.  However, it illustrates the improvement due to MLEs later:
%STL1Real.makeFitPlot

%% Compute the MLE
STL1pars = STL1Real.maximizeLikelihood(STL1pars, fitOptions);

% Make plots of the model parameter fits from the MLEs
STL1Real.makeFitPlot

%% Compute the FIM
Model_fimResults = ModelReal.computeFIM([],'log');
STL1_fimResults = STL1Real.computeFIM([],'log'); 

% Generate a count of measured cells (i.e., in the case of real data,  
% the number of cells measured in the experiment)
Model_cellCounts = 10*ones(size(ModelReal.tSpan));
STL1_cellCounts = 10*ones(size(STL1Real.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[Model_fimTotal,Model_mleCovEstimate,Model_fimMetrics] = ...
     ModelReal.evaluateExperiment(Model_fimResults,Model_cellCounts)

[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
     STL1Real.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

% Plot the FIMs
fig1 = figure(1);clf; set(fig1,'Name',...
    'Fim-Predicted Uncertainty Ellipses');
ModelReal.plotMHResults([],Model_fimTotal,'log',[],fig1)
legend('FIM')

fig2 = figure(2);clf; set(fig2,'Name',...
    'Fim-Predicted Uncertainty Ellipses');
STL1Real.plotMHResults([],STL1_fimTotal,'log',[],fig2)
legend('FIM')