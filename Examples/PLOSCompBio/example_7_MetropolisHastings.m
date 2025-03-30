%% example_7_MetropolisHastings
% Example script to demonstrate how to use Metropolis-Hastings to sample
% uncertainty
close all
addpath(genpath('../src'));

%% Preliminaries
% Load our models described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
example_6_LoadingandFittingData_MLE
Model.summarizeModel
STL1Model.summarizeModel