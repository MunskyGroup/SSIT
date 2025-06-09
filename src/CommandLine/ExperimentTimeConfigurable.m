classdef ExperimentTimeConfigurable < AbstractExperimentConfigurable
    properties      
        Values (1, :) double {mustBeNonempty, mustBeNonnegative} = [0]
    end
    methods
        function disp(obj)           
            disp(['Time = ' num2str(obj.Values)])            
        end
        function value = getSingleValue(obj)
            if length(obj.Values) > 1
                error('Called getSingleValue on configurable ' + obj.disp)
            else
                value = obj.Values(1);
            end
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
    methods (Static)
        function varName = getVarName()
            varName = 'Time';
        end
        function varType = getVarType()
            varType = 'double';
        end
    end
end