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

    methods (Access = protected)
        function isEqual = equalsSibling(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentInputConfigurable
                obj2 (1, 1) MBExperimentInputConfigurable
            end

            % The Values setter guarantees that all values will be unique
            % and in sorted order, so element-wise comparison is fair:

            isEqual = strcmp(obj1.InputName, obj2.InputName) && ...
                length(obj1.Values) == length(obj2.Values) && ...
                all(strcmp(obj1.Values, obj2.Values));            
        end % equalsSibling
    end % Protected instance methods

    methods
        function model = applyToModel(obj, model)
            %applyToModel takes a model and configures it according to this
            %configurable, provided that it has only a single value, since
            %SSIT does not allow a single input expression to take multiple
            %values.           
            
            arguments
                obj (1, 1) MBExperimentInputConfigurable ...
                    {mustHaveSingleValue}                
                model (1, 1) MBExperimentModel
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
            fprintf("%s = %s\n", ...
                obj.InputName, join(string(obj.Values), ", "))           
        end

        function rows = findSubsetOfData(obj, data)
            %findSubsetOfData takes a dataset in the form of a MATLAB table
            %and finds its subset that matches this configurable; in other
            %words, it will return a logical row vector containing true
            %exactly where the data table contains a row for which the
            %value of the column corresponding to the InputName equals one
            %of the Values.

            arguments
                obj (1, 1) MBExperimentInputConfigurable
                data (1, 1) table
            end

            % We use strcmp because the values are strings and potentially
            % case-sensitive.

            rowsToKeep = strcmp(data.(obj.InputName), obj.Values(1));

            for valIdx = 2:length(obj.Values)
                % If there are multiple values, then each row that matches
                % any of the values should be kept. This is equivalent to
                % logical disjunction ("or").

                curRowsToKeep = ...
                    strcmp(data.(obj.InputName), obj.Values(valIdx));
                rowsToKeep = or(rowsToKeep, curRowsToKeep);
            end

            % Return the identified subset of the data rows.

            rows = rowsToKeep;
        end % findSubsetOfData

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