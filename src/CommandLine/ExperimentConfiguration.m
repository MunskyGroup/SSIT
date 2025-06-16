classdef ExperimentConfiguration
    properties
        NumberOfObservations (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
        Configurables (1, :) AbstractExperimentConfigurable = createArray(1,1,"ExperimentTimeConfigurable")
    end
    methods
        function disp(obj)
            for configIdx = 1:length(obj.Configurables)
                disp(obj.Configurables(configIdx))
            end
            disp([num2str(obj.NumberOfObservations), ' Observation(s)'])            
        end
    
        function model = applyToModel(obj, model)
            arguments
                obj
                model (1, 1) DiscoverableModel
            end

            for configIdx = 1:length(obj.Configurables)
                model = obj.Configurables(configIdx).applyToModel(model);
            end
        end
    end
end