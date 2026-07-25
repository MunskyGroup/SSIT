classdef MBExperimentDesignStrategyInfo
    %MBExperimentDesignStrategyInfo contains information needed by a
    %specific strategy for the latter to design an experiment. This
    %ancestor contains information needed by all strategies; descendants
    %can provide additional information needed by specific strategies.

    properties
        NumberOfObservationsPerExperiment (1, 1) double ...
            {mustBeInteger, mustBePositive} = 10
        ObservationQuantum (1, 1) double ...
            {mustBeInteger, mustBePositive} = 1
    end

    methods
        function mustBeValid(obj)
            arguments
                obj (1, 1) MBExperimentDesignStrategyInfo
            end

            if any(mod(obj.NumberOfObservationsPerExperiment, ...
                    obj.ObservationQuantum))
                error("Attempt to design experiment whose " + ...
                    "observation total (" + ...
                    num2str(obj.NumberOfObservationsPerExperiment) + ...
                    ") is not an even multiple of the quantum (" + ...
                    num2str(obj.ObservationQuantum) + ")")
            end
        end
    end
end