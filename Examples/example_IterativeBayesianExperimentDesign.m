close all
clear all
addpath('../CommandLine')

%% Define Model
ModelTrue = SSIT;
ModelTrue.species = {'rna'};
ModelTrue.initialCondition = 0;
ModelTrue.propensityFunctions = {'kr';'gr*rna'};
ModelTrue.stoichiometry = [1,-1];
ModelTrue.fspOptions.fspTol = 1e-5;
ModelTrue.parameters = ({'kr',10;'gr',0.3});
ModelTrue = ModelTrue.formPropensitiesGeneral('Poiss',true);

%% Assign 'true' parameters
trueParameters = [10;0.3];
ModelTrue.parameters = ({'kr',trueParameters(1);'gr',trueParameters(2)});
ModelTrue = ModelTrue.formPropensitiesGeneral('Poiss',true);
[ModelSolution,Model.fspOptions.bounds] = ModelTrue.solve;

%% Possible Time points for Experiments.
nT = 21;
ModelTrue.tSpan = linspace(0,10,nT);

%% Set up for first experiment
numCellsPerExperiment = 15;  % Total number of cells in each experiment
nextExperiment = zeros(1,nT);
nextExperiment([1,11,21]) = round(numCellsPerExperiment/3);

%% Loop over subsequent experiments
Nexp = 5;
for iExpt = 1:5
    % Generate "Fake" data
    dataFile = ['FakeExperiment_',num2str(iExpt),'.csv'];
    ModelTrue.sampleDataFromFSP(testCase.PoissSolution,'testData.csv')



    % 


end


