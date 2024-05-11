% Gene Expression System Controller Class Definition
classdef (Abstract) GeneExprSysController < handle
    properties (SetAccess = private)
        CompletedRounds (1,1) integer
        EstimatedModel (1,1) GeneExprSysModel
        InputLibrary (1,:) cell
        MaxRounds (1,1) integer {mustBePositive}
        SaveName (1,1) string {mustBeNonempty}
    end

    methods
        function c = GeneExprSysController(model, rounds, inputs, saveName)
            c.CompletedRounds = 0;
            c.EstimatedModel = model;
            if isempty(inputs)
                warning('Controller Inputs Cannot Be Empty')
                return
            end
            c.InputLibrary = inputs;
            c.MaxRounds = rounds;
            c.SaveName = lower([saveName,'.mat']);
            % Exit if the savefile already exists.
            if exist(c.SaveName,'file')
                warning('Save File Already Exists -- Skipping Control')
                return
            end
        end
    end

    methods(Abstract, Access = ?GeneExprSysController)
        experiment = SelectNextExperiment(controller)
        updateModel(controller, data) 
    end
end