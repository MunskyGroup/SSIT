classdef MBExperimentUniformDesignStrategy < ...
        MBAbstractExperimentDesignStrategy

    methods (Access = protected)
        function [design, cellVecForOptimalFIMCalculation] = ...
                apportionObservationsInternal(obj, design)
            arguments
                obj (1, 1) MBExperimentUniformDesignStrategy
                design (1, 1) MBExperimentDesign
            end

            % We will ignore the info's ObservationQuantum here; that
            % property is intended to simplify more complicated design
            % strategies. Instead, we will simply give each configuration
            % the largest number of observations possible without exceeding
            % the total limit:

            observationsPerConfiguration = ...
                floor(obj.Info.NumberOfObservationsPerExperiment / ...
                design.NumberOfConfigurations);
            configurations = design.Configurations;
            for configIdx = 1:length(configurations)
                design.setObservationsForConfiguration(...
                    configurations(configIdx), ...
                    observationsPerConfiguration);
            end

            [Obs, ~, ~] = design.getAsObservationMatrix();
            Obs = reshape(Obs', [1, numel(Obs)]);
            cellVecForOptimalFIMCalculation = Obs';
        end
    end
end