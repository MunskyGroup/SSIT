% Run this script when paused on line 427 in 
% iterativeExperimentRunnerMultiConditions.m
% Note that the runner needs to be temporarily edited to remove these
% class constraints: MBExperimentModel, MetropolisHastingsAlgorithmOptions 

runner = MetropolisHastingsAlgorithmRunner();
runner.FIMTrue = fimTrue;
runner.FitOptions = fitOptions;
runner.GuessedModelsWithCombinedTimes = ModelGuess;

runner.IDString = datFileName;
runner.IDNumber = ind;
MHFitOptions.Progress = false;
MHFitOptions.SaveFile = ['TMPMHChain_',datFileName,'_',num2str(ind),'.mat'];
runner.MHOptions = MHFitOptions;
runner.ModelsHaveData = hasData;
runner.NumberOfSamplesForProduction = min(nSamplesMH, 200);
runner.ObservationsMatrix = nTotalCells;
runner.ParameterGuesses = newPars;
runner.StateSpaces = stateSpaces;
runner.TrueModelsWithCombinedTimes = ModelTrue;
runner.UsingSimulatedData = true; % default

[results, stateSpaces] = runner.run();
results