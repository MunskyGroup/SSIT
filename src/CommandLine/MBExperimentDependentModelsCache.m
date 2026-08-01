classdef MBExperimentDependentModelsCache
    %MBExperimentDependentModelsCache holds a set of models and related
    %properties, with the intention that an MBExperimentDesigner can use it
    %to remember a set of models created earlier for some particular
    %purpose (e.g., to fit models to data) and quickly determine whether
    %the models are up-to-date, saving computational time by avoiding
    %unnecessary recalculations.
    %
    %No user should ever interact with this class directly; thus, all
    %property access is restricted to MBExperimentDesigner.

    properties (Access = ?MBExperimentDesigner)
        Configs (1, :) MBExperimentConfiguration
        LastPerformedExperimentRound (1, 1) uint64 = 0
        Models (1, :) MBExperimentModel
        ModelsHaveData (1, :) logical
        ModelStateSpaces (1, :) cell 
        SourceModel (1, 1) MBExperimentModel
        UseEmpiricalData (1, 1) logical = false
    end

    methods
        function obj = untitled(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end