classdef DesignedExperiment
    properties
        % The BaseModel should be assigned by the object owning the
        % experiment. It reflects the model as it should be used EXCEPT FOR
        % any modification(s) that the configuration might make.        
        BaseModel (1, 1) DiscoverableModel
        Configuration (1, 1) ExperimentConfiguration       
    end

    properties (Dependent)
        % The Model should be queried by any object referencing the
        % experiment, e.g. data-fitting and simulation routines. It 
        % includes any modification(s) that the configuration might make.
        % Being a dependent property, it is recalculated upon demand to
        % reflect the latest values of both the base model and the
        % configuration.
        Model (1, 1) DiscoverableModel
    end

    methods
        function model = get.Model(obj)
            model = obj.Configuration.applyToModel(obj.BaseModel);
        end
    end
end