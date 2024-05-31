% Gene Expression System Controller Class Definition
classdef (Abstract) GeneExprSysController < handle
    properties (SetAccess = private)
        MaxRounds (1,1) integer {mustBePositive}
        SaveName (1,1) string {mustBeNonempty}
    end
    properties (Access = protected)
        AllDataSoFar (1,:) cell
        CompletedRounds (1,1) integer
        NumCellsPerExperiment (1,1) integer {mustBePositive}
        InputLibrary (1,:) cell
        EstimatedModel (1,1) GeneExprSysModel
    end

    methods
        function c = GeneExprSysController(...
                model, rounds, inputs, cellsPerExperiment, saveName)
            c.CompletedRounds = 0;
            c.EstimatedModel = model;
            if isempty(inputs)
                warning('Controller Inputs Cannot Be Empty')
                return
            end
            c.AllDataSoFar = cell(1,length(inputs));
            c.InputLibrary = inputs;
            c.NumCellsPerExperiment = cellsPerExperiment;
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