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
        % evaluateExperiment expects a column vector of FIMs, so it is most
        % efficient for us to convert it here, since we will pull from the
        % cache whenever possible.
        FIMs (:, 1) cell
        LastPerformedExperimentRound (1, 1) uint64 = 0
        Models (1, :) MBExperimentModel
        ModelsHaveData (1, :) logical
        ModelStateSpaces (1, :) cell 
        SourceModel (1, 1) MBExperimentModel
        UseEmpiricalData (1, 1) logical = false
    end
end