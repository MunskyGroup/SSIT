classdef ExperimentInputConfigurable < AbstractExperimentConfigurable
    properties
        InputName (1, 1) string
        Values (1, :) double {mustBeNonempty} = [0]
    end
    methods
        function disp(obj)
            fprintf('%s = %s\n', obj.InputName, num2str(obj.Values))           
        end
        function value = getSingleValue(obj)
            if length(obj.Values) > 1
                error('Called getSingleValue on configurable ' + obj.disp)
            else
                value = obj.Values(1);
            end
        end
        function varName = getVarName(obj)
            varName = obj.InputName;
        end
        function configurables = multiply(obj)
            configurables = createArray(1, length(obj.Values), ...
                "ExperimentInputConfigurable");
            for valuesIdx = 1:length(obj.Values)
                configurables(valuesIdx).InputName = obj.InputName;
                configurables(valuesIdx).Values = obj.Values(valuesIdx);               
            end
        end
        function possibilities = numberOfValues(obj)
            possibilities = length(obj.Values);
        end
    end
    methods (Static)
        function varType = getVarType()
            varType = 'double';
        end
    end
end