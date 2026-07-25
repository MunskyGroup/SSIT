classdef MBExperimentRound
    %MBEXPERIMENTROUND This class holds the summary of a single round of
    %experiment design.

    properties
        CovarianceMH
        CovarianceLogMH        
        CovarianceFIM_Prediction
        CumulativeNumbersOfObservations (1, :) uint64 % Before the round
        FIMCurrentExpt
        FIMCurrentExptTrue
        FIMOptNextExpt cell
        FIMPredNextExpt
        FIMResults
        MHResults (1, 1) MetropolisHastingsAlgorithmResults 
        NextExperimentDesign (1, 1) MBExperimentDesign
        ParametersFound
    end

    methods
    end
end