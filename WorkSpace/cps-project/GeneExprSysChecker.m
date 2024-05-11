% Gene Expression System Checker Class Definition
classdef GeneExprSysChecker < matlab.mixin.Heterogeneous
    % Inheriting from matlab.mixin.Heterogeneous means you can make arrays
    % of objects from the subclasses, with different subclasses in the same
    % array.

    properties (SetAccess = private)
        BaseModel (1,1) GeneExprSysModel
        Controller (1,1) GeneExprSysController
        Invariants (:,1) Invariant
        Plant (1,1) GeneExprSysPlant
    end

    methods
        function c = GeneExprSysChecker(model, controller, invariants, plant)
            c.BaseModel = model;
            c.Controller = controller;
            c.Invariants = invariants;
            c.Plant = plant;
        end
    end

    methods(Sealed)
        % Sealed methods cannot be redefined in subclasses
        function valid = check_system(checker)
            valid = true;

            % Perform experiments while possible
            while checker.Controller.CanExperiment()
                nextExperiment = checker.Controller.SelectNextExperiment();
                latestData = checker.Plant.RunNextExperiment(nextExperiment);
                checker.Controller.UpdateModel(latestData);
            end

            % Return true if and only if all invariants are satisfied
            for i = 1:length(checker.Invariants)
                valid = checker.Invariants(i).predicate();
                if ~valid
                    break
                end
            end
        end
    end
end