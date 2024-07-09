function [saveExpName] = iterativeExperimentRunnerMultiConditions(example,data,sampleType,nExptRounds,...
    rngSeed,ind,initialExperiment,...
    nFIMsamples,truePars,saveFileName,initialParameterGuess,...
    testing)

% Default MLE Fitting Options
maxFitIter = 200;
nFitRounds = 3;

% Maximum number of cells is infinite unless specified otherwise.
maxAvailable = []; 

ModelSolution = cell(1,nInputs);

% Randomize the initial parameter set.
newPars = initialParameterGuess;

%% Loop over subsequent experiments
% Keep track of accumulated data in each experiment design
  
nTotalCells = zeros(nInputs,nT); 

for iExpt = 1:nExptRounds
    nTotalCells = nTotalCells + nextExperiment;
 

    

    
    
    

    %% MH Tuning then running.

    

end % for iExpt