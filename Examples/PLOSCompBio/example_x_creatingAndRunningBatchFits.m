%% Example -- Running SSIT Pipelines Outside of Matlab
close all
clear 
addpath(genpath('../../src'));
%% Define Model/Data Combination
% Specify data set to fit.
DataSettings = {'data/STL1.csv',{'mRNA','rna'}};

% Create model from preset, associate with data.
Model = SSIT('BirthDeath',[],DataSettings);
Model = Model.formPropensitiesGeneral('BirthDeath');

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

%% Run on cluster
% Note that this is set up for a specific cluster that uses sbatch. To run
% on your own cluster, you will need to change the argument
% 'clusterPrefix' below.
% Generate the command so you can run it from terminal.
% cmd = SSIT.generateCommandLinePipeline(saveFile,'Model',[],Pipeline,pipelineArgs,saveFile,'logFile.txt',0,1)

% Run it directly.
% cmd = SSIT.generateCommandLinePipeline(saveFile,'Model',[],Pipeline,pipelineArgs,saveFile,'logFile.txt',1,1)


