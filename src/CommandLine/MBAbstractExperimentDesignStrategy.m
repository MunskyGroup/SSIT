classdef MBAbstractExperimentDesignStrategy
    %MBExperimentDesignStrategy determines how to apportion available
    %observations across candidate experiment configurations.

    properties
        Info (1, 1) MBExperimentDesignStrategyInfo {mustBeValid(Info)}
    end

    methods (Abstract, Access = protected)
        [design, cellVecForOptimalFIMCalculation] = ...
            apportionObservationsInternal(obj, design)                
    end % Abstract protected methods

    methods (Access = protected)
        function obj = incorporateRoundDetailsInternal(obj, round)
            arguments
                obj (1, 1) MBAbstractExperimentDesignStrategy
                round (1, 1) MBExperimentRound
            end

            % Do nothing here - descendants can override.
        end
    end % Protected methods

    methods
        function [design, cellVecForOptimalFIMCalculation] = ...
                apportionObservations(obj, design)
            arguments (Input)
                obj (1, 1) MBAbstractExperimentDesignStrategy
                design (1, 1) MBExperimentDesign
            end

            [design, cellVecForOptimalFIMCalculation] = ...
                apportionObservationsInternal(obj, design);
            mustNotExceedMaximumObservations(design, obj);
        end

        function obj = incorporateRoundDetails(obj, round)
            arguments
                obj (1, 1) MBAbstractExperimentDesignStrategy
                round (1, 1) MBExperimentRound
            end

            obj = incorporateRoundDetailsInternal(obj, round);
        end
    end % Public methods
end

function mustNotExceedMaximumObservations(design, strategy)
arguments
    design (1, 1) MBExperimentDesign
    strategy (1, 1) MBAbstractExperimentDesignStrategy
end

if design.NumberOfObservations > ...
        strategy.Info.NumberOfObservationsPerExperiment
    error("Too many observations (" + ...
        num2str(design.NumberOfObservations) + ...
        ") were apportioned to an experiment with a maximum of " + ...
        num2str(strategy.Info.NumberOfObservationsPerExperiment))
end

end