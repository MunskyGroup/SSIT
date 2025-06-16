classdef ExperimentInputConfigurable < AbstractExperimentConfigurable
    properties
        InputName (1, 1) string
        Values (1, :) double {mustBeNonempty} = [0]
    end
    
    methods
        function model = applyToModel(obj, model)
            arguments
                obj                 
                model (1, 1) DiscoverableModel 
            end

            % We can only apply an input configurable to a model when it
            % has a single value; calling getSingleValue here does
            % double-duty, since it will cause an error if there are
            % multiple Values.

            value = obj.getSingleValue();

            % We apply the input configurable to the model by concatenating
            % the appropriate pair to its inputExpressions array. Note that
            % by convention, the model uses an 'I' prefix for input
            % variables, so we prepend that here. Also, we express all
            % numerical values as strings within input expressions, so that
            % conversion is required here.

            model.inputExpressions = [model.inputExpressions; ...
                {append('I', obj.InputName), num2str(value)}];
        end % applyToModel

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