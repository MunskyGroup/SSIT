%% example_16_SSITPipelinesAndClusterComputing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, data loaded in 
% example_8_LoadingandFittingData_DataLoading, and MLE computed in
% example_9_LoadingandFittingData_MLE
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE

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
cmd = SSIT.generateCommandLinePipeline(saveFile,'Model',[],Pipeline,pipelineArgs,saveFile,'logFile.txt',0,1)

% Run it directly.
cmd = SSIT.generateCommandLinePipeline(saveFile,'Model',[],Pipeline,pipelineArgs,saveFile,'logFile.txt',1,1)


%% Running a batch job for cross validation
%%      Example 6 -- using multimodel to explore batch variations
% In this example, we will use the batch job capabilities to do
% cross-validation where we fit different replicas at the same time.

% Let's first create a new model.
Model1 = SSIT;
Model1.species = {'onGene';'rna'};
Model1.initialCondition = [0;0];
Model1.propensityFunctions = {'kon*IGR*(2-onGene)';'koff*onGene';'kr*onGene';'gr*rna'};
Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model1.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
Model1.parameters = ({'koff',0.14;'kon',0.14;'kr',10;'gr',0.01;...
    'a1',0.4;'r1',0.04;'r2',0.1});
Model1.fspOptions.initApproxSS = true;
Model1.fittingOptions.modelVarsToFit = 1:7;
% We generate functions for model propensities
Model1 = Model1.formPropensitiesGeneral('Model1FSP');

% Specify datafile name and species linking rules
DataFileName = '../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv';
LinkedSpecies = {'rna','RNA_nuc'};

% In this case, let's suppose that we only wish to fit the data at times
% before 75 minutes.  We will set the global conditions as:
ConditionsGlobal = {[],[],'TAB.time<=75'};

% We want to split up the replicas to be separate.
ConditionsReplicas = {'TAB.Rep_num==1';...
    'TAB.Rep_num==2'};

modelLibrary = 'crossValidationModelsDusp1';

SSIT.runCrossValidation(Model1,DataFileName, ...
    LinkedSpecies,ConditionsGlobal,ConditionsReplicas, ...
    modelLibrary)