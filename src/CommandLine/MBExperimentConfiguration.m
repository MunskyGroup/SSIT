classdef MBExperimentConfiguration
    %MBExperimentConfiguration encapsulates a configuration, i.e. an
    %assignment of a single value for all "configurables": the non-time
    %configurations (input expressions and designed parameters of a model)
    %and the measurement time(s).
    
    properties
        NonTimeConfiguration (1, 1) MBExperimentNonTimeConfiguration
        TimeConfigurable (1, 1) MBExperimentTimeConfigurable
    end

    methods
        function model = applyToModel(obj, model)
            %applyToModel takes an SSIT model and configures it according
            %to these configurables.
            arguments
                obj                 
                model (1, 1) SSIT
            end

            nonTimeConfigurables = ...
                obj.NonTimeConfiguration.NonTimeConfigurables;
            for nonTimeConfigIdx = 1:length(nonTimeConfigurables)
                curConfig = nonTimeConfigurables(nonTimeConfigIdx);
                model = curConfig.applyToModel(model);
            end

            model = obj.TimeConfigurable.applyToModel(model);
        end % applyToModel
    end
end