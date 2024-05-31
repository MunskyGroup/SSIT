% Gene Expression System Plant Class Definition
classdef (Abstract) GeneExprSysPlant < handle
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