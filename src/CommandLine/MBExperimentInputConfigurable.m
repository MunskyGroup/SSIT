classdef MBExperimentInputConfigurable < ...
        MBAbstractExperimentNonTimeConfigurable
    %MBExperimentInputConfigurable represents an individual input
    %configurable, i.e., an input expression with a specific name and one
    %or more possible values.
    %
    %   Example:
    %       InputName = "IDex"
    %       Values = ["0" "10" "20" "30*(t<30)+50*(t>=30)"]
    %   In this example, the "IDex" input to the model can take one of four
    %   possible values. Observe that the value can be time-varying and can
    %   include arithmetical, logical, and relational operators.

    properties
        InputName (1, 1) string {mustBeNonempty} = "IInput"
        Values (1, :) string {mustBeNonempty} = "0"
    end

    properties (Dependent)
        FilenameString
    end

    methods
        function model = applyToModel(obj, model)
            %applyToModel takes an SSIT model and configures it according
            %to this configurable, provided that it has only a single
            %value, since SSIT does not allow a single input expression to
            %take multiple values.           
            
            arguments
                obj (1, 1) MBExperimentInputConfigurable ...
                    {mustHaveSingleValue}                
                model (1, 1) SSIT
            end

            value = obj.getSingleValue();

            % We apply the input configurable to the model by concatenating
            % the appropriate pair to its inputExpressions array; because
            % the latter historically uses char arrays for both "key" and
            % "value," we perform that conversion here.

            model.inputExpressions = [model.inputExpressions; ...
                {char(obj.InputName), char(value)}];
        end % applyToModel

        function disp(obj)
            fprintf('%s = %s\n', obj.InputName, obj.Values)           
        end

        function filenameString = get.FilenameString(obj)
            %FilenameString represents the configurable as a string
            %suitable for use in a filename.

            filenameString = join([obj.InputName obj.Values], "_");
        end

        function value = getSingleValue(obj)
            arguments
                obj (1, 1) MBExperimentInputConfigurable ...
                    {mustHaveSingleValue}
            end
            
            value = obj.Values(1);        
        end

        function varName = getVarName(obj)
            varName = obj.InputName;
        end

        function result = lt(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentInputConfigurable
                obj2 (1, 1) MBExperimentInputConfigurable
            end

            if obj1.InputName < obj2.InputName
                result = true;
            elseif obj1.InputName > obj2.InputName
                result = false;
            else
                % The input names are the same. The smaller one will be the
                % one that has fewer values. If each object has the same
                % number of values, iterate over each pair of elements to
                % determine the smaller.
                if obj1.numberOfValues() < obj2.numberOfValues()
                    result = true;
                elseif obj1.numberOfValues() > obj2.numberOfValues()
                    result = false;
                else
                    result = false;
                    for valIdx = 1:obj1.numberOfValues()
                        if obj1.Values(valIdx) < obj2.Values(valIdx)
                            result = true;
                            break;
                        elseif obj1.Values(valIdx) > obj2.Values(valIdx)
                            result = false;
                            break;
                        end
                    end
                end
            end
        end % lt

        function configurables = multiply(obj)
            %multiply creates an array of input configurables each
            %containing a distinct single value in the original input
            %configurable.
            configurables = createArray(1, length(obj.Values), ...
                "MBExperimentInputConfigurable");
            for valuesIdx = 1:length(obj.Values)
                configurables(valuesIdx).InputName = obj.InputName;
                configurables(valuesIdx).Values = obj.Values(valuesIdx);               
            end
        end

        function possibilities = numberOfValues(obj)            
            possibilities = length(obj.Values);
        end

        function obj = set.Values(obj, val)
            % For simplicity and clarity, we will remove duplicate
            % expressions and sort them in ascending order.

            obj.Values = unique(sort(val));            
        end
    end % Public instance methods

    methods (Static)
        function varType = getVarType()
            varType = "string";
        end
    end
end

function mustHaveSingleValue(configurable)
arguments
    configurable (1, 1) MBExperimentInputConfigurable
end

if configurable.numberOfValues > 1
    error("An input configurable " + configurable.FilenameString + ...
        " contained multiple values where exactly one is required.");
end

end