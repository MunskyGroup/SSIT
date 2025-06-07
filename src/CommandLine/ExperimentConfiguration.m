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
    end
end