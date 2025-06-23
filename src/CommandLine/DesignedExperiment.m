classdef DesignedExperiment
    properties
        Configuration (1, 1) ExperimentConfiguration       
        % The UnderlyingModel should be assigned by the object owning the
        % experiment. It reflects the model as it should be used EXCEPT FOR
        % any modification(s) that the configuration might make.        
        UnderlyingModel (1, 1) DiscoverableModel        
    end

    properties (Dependent)
        % The Model should be queried by any object referencing the
        % experiment, e.g. data-fitting and simulation routines. It 
        % includes any modification(s) that the configuration might make.
        % Being a dependent property, it is recalculated upon demand to
        % reflect the latest values of both the underlying model and the
        % configuration.
        Model (1, 1) DiscoverableModel
    end

    methods
        function model = get.Model(obj)
            model = obj.Configuration.applyToModel(obj.UnderlyingModel);
        end
    end
end