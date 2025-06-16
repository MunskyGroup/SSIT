classdef ExperimentTimeConfigurable < AbstractExperimentConfigurable
    properties      
        Values (1, :) double {mustBeNonempty, mustBeNonnegative} = [0]
    end
    
    methods
        function model = applyToModel(obj, model)
            arguments
                obj                 
                model (1, 1) DiscoverableModel 
            end

            % Although in principle we could apply a multi-valued time 
            % configurable to a model, for now we will not do so; calling 
            % getSingleValue here does double-duty, since it will cause an 
            % error if there are multiple Values.

            value = obj.getSingleValue();

            % We apply the time configurable to the model by setting the
            % tSpan property of the latter accordingly.

            model.tSpan = value;
        end % applyToModel

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
    end % Public instance methods

    methods (Static)
        function varName = getVarName()
            varName = 'Time';
        end

        function varType = getVarType()
            varType = 'double';
        end
    end
end