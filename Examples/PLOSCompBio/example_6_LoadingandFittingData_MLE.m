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
STL1Model.summarizeModel

% Set fitOptions, with the maximum allowable number of iterations to fit
fitOptions = optimset('Display','iter','MaxIter',500);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
STL1pars = cell2mat(STL1Model.parameters(1:8,2));

%% Load experimental data, matching the species name of the model ('mRNA') 
% to the appropriate column in the data file ('rna')
STL1Real = STL1Model.loadData('data/STL1.csv',{'mRNA','rna'});

%% Compute the MLE
STL1pars = STL1Real.maximizeLikelihood(STL1pars, fitOptions);