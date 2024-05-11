% Gene Expression System Plant Class Definition
classdef (Abstract) GeneExprSysPlant < matlab.mixin.Heterogeneous
    % Inheriting from matlab.mixin.Heterogeneous means you can make arrays
    % of objects from the subclasses, with different subclasses in the same
    % array.

    properties (SetAccess = private)
        InputLibrary (1,:) cell
        TrueModel (1,1) GeneExprSysModel
    end

    methods
        function p = GeneExprSysPlant(model, inputs)
            if isempty(inputs)
                warning('Plant Inputs Cannot Be Empty')
                return
            end
            p.InputLibrary = inputs;
            p.TrueModel = model;
        end
    end

    methods(Abstract)
        data = RunNextExperiment(plant, experiment)
    end
end