classdef ExperimentTimeConfigurable < AbstractExperimentConfigurable
    properties      
        Values (1, :) double {mustBeNonempty, mustBeNonnegative} = [0]
    end
    methods
        function disp(obj)           
            disp(['Time = ' num2str(obj.Values)])            
        end
        function configurables = multiply(obj)
            configurables = createArray(1, length(obj.Values), ...
                "ExperimentTimeConfigurable");
            for valuesIdx = 1:length(obj.Values)             
                configurables(valuesIdx).Values = obj.Values(valuesIdx);               
            end
        end
        function possibilities = numberOfValues(obj)
            possibilities = length(obj.Values);
        end
    end
end