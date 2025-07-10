%% Example -- Running SSIT Pipelines Outside of Matlab
close all
clear 
addpath(genpath('../../src'));
%% Define Model/Data Combination
% Specify data set to fit.
DataSettings = {'data/STL1.csv',{'mRNA','rna'}};

% Create model from preset, associate with data.
Model = SSIT('BirthDeath',[],DataSettings);

% Save model for later use.
saveFile = 'generatedModel.mat';
save(saveFile,"Model")

%% Call Pipeline to Fit Model
% Specify pipeline to apply to model and arguments
% ("../src/exampleData/examplePipelines/fittingPipelineExample.m") 
Pipeline = 'fittingPipelineExample';
pipelineArgs.maxIter = 1000;
pipelineArgs.display = 'iter';
pipelineArgs.makePlot = false;

SSIT(saveFile,'Model',[],Pipeline,pipelineArgs,saveFile);

%% Get command to run job in background.
% Here is a way to generate the command so you can run it from terminal.
cmd = SSIT.generateCommandLinePipeline(saveFile,'Model',[],Pipeline,pipelineArgs,saveFile,'logFile.txt')

%% Run that command from matlab.
system(cmd);