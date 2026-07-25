classdef MBExperimentRandomDesignStrategy < ...
        MBAbstractExperimentDesignStrategy

    methods (Access = protected)
        function [design, cellVecForOptimalFIMCalculation] = ...
                apportionObservationsInternal(obj, design)
            arguments
                obj (1, 1) MBExperimentRandomDesignStrategy
                design (1, 1) MBExperimentDesign
            end

            quanta = obj.Info.NumberOfObservationsPerExperiment / ...
                obj.Info.ObservationQuantum;

            apportionToConfigs = randi(design.NumberOfConfigurations, ...
                1, quanta);

            for configIdx = apportionToConfigs
                curConfig = design.Configurations(configIdx);
                curNumberOfObservations = ...
                    design.getObservationsForConfiguration(curConfig);
                design.setObservationsForConfiguration(curConfig, ...
                    curNumberOfObservations + obj.Info.ObservationQuantum);
            end

            [Obs, ~, ~] = design.getAsObservationMatrix();
            Obs = reshape(Obs', [1, numel(Obs)]);
            cellVecForOptimalFIMCalculation = Obs';            
        end % apportionObservationsInternal
    end % Protected methods

    methods
        function obj = MBExperimentRandomDesignStrategy(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Info = MBExperimentDesignStrategyInfo;
            obj.Property1 = inputArg1 + inputArg2;
        end
    end % Public methods
end