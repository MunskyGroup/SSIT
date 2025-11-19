%% example_16_PipelinesAndClusterComputing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.7: Pipelines and Cluster Computing Scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, data loaded in 
% example_8_LoadingandFittingData_DataLoading, and MLE computed in
% example_9_LoadingandFittingData_MLE
%clear
%close all
addpath(genpath('../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE
% example_10_LoadingandFittingData_MH

%% Load model fitted using Metropolis-Hastings:
% load('example_10_LoadingandFittingData_MH.mat')

% View summary of 4-state STL1 model:
STL1_4state_MH.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running SSIT Pipelines Outside of Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of our 4-state STL1 model:
STL1_4state_pipe = STL1_4state_MH;

% Save model for later use:
saveFile = 'generatedModel.mat';
save(saveFile,"STL1_4state_pipe")

%% Call Pipeline to Fit Model
% Specify pipeline to apply to model and arguments
% ("../src/exampleData/examplePipelines/fittingPipelineExample.m") 
Pipeline = 'fittingPipelineExample';
pipelineArgs.maxIter = 1000;
pipelineArgs.display = 'iter';
pipelineArgs.makePlot = false;

SSIT(saveFile,'STL1_4state_pipe',[],Pipeline,pipelineArgs,saveFile);

%% Get command to run job in background.
% Here is a way to generate the command so you can run it from terminal.
cmd = SSIT.generateCommandLinePipeline(saveFile,'STL1_4state_pipe',[],Pipeline,pipelineArgs,saveFile,'logFile.txt')

%% Run that command from matlab.
system(cmd);

%% Run on cluster
% Note that this is set up for a specific cluster that uses sbatch. To run
% on your own cluster, change the argument 'clusterPrefix' 
% Generate the command so you can run it from terminal.
cmd = SSIT.generateCommandLinePipeline(saveFile,'STL1_4state_pipe',[],Pipeline,pipelineArgs,saveFile,'logFile.txt',0,1)

% Run it directly.
cmd = SSIT.generateCommandLinePipeline(saveFile,'STL1_4state_pipe',[],Pipeline,pipelineArgs,saveFile,'logFile.txt',1,1)