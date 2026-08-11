classdef MBExperimentRound
    %MBEXPERIMENTROUND This class holds the summary of a single round of
    %experiment design.

    properties
        CovarianceMH
        CovarianceLogMH        
        CovarianceFIM_Prediction
        % The cumulative experiment design refers to all PERFORMED rounds 
        % preceding this DESIGNED one.
        CumulativeExperimentDesign (1, 1) MBExperimentDesign
        FIMCurrentExpt
        FIMCurrentExptTrue
        FIMOptNextExpt cell
        FIMPredNextExpt
        FIMResults
        FIMTrue
        MHResults (1, 1) MetropolisHastingsAlgorithmResults 
        NextExperimentDesign (1, 1) MBExperimentDesign
        ParametersFound
    end

    methods
    end
end