classdef MBExperimentTimeConfigurable < MBAbstractExperimentConfigurable
    %MBExperimentTimeConfigurable represents an individual time
    %configurable, i.e., one or more possible system
    %observation/measurement times.
    %
    %   Example:    
    %       Values = [0 10.5 21 31.5 50 90]
    %   In this example, the system can be observed and measurements made
    %   at any of six distinct times. Observe that the values need not be
    %   integral nor evenly spaced; however, they are restricted to
    %   non-negative values.

    properties      
        Values (1, :) double {mustBeNonempty, mustBeNonnegative} = [0]
    end

    properties (Dependent)
        FilenameString
    end

    methods (Access = protected)
        function isEqual = equalsSibling(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentTimeConfigurable
                obj2 (1, 1) MBExperimentTimeConfigurable
            end

            % The Values setter guarantees that all values will be unique
            % and in sorted order, so element-wise comparison is fair:

            isEqual = length(obj1.Values) == length(obj2.Values) && ...
                all(obj1.Values == obj2.Values);            
        end % equalsSibling
    end % Protected instance methods

    methods
        function model = applyToModel(obj, model)
            %applyToModel takes a model and configures it according to this
            %configurable. SSIT allows a time span to have multiple values,
            %and the setter of the Values property already ensures that all
            %times are unique and in ascending order.
            
            arguments
                obj (1, 1) MBExperimentTimeConfigurable
                model (1, 1) MBExperimentModel
            end

            % We apply the time configurable to the model by setting the
            % tSpan property of the latter accordingly.

            model.tSpan = obj.Values;
        end % applyToModel

        function disp(obj)           
            disp("Time = " + join(string(obj.Values), ", "))            
        end

        function rows = findSubsetOfData(obj, data)
            %findSubsetOfData takes a dataset in the form of a MATLAB table
            %and finds its subset that matches this configurable; in other
            %words, it will return a logical row vector containing true
            %exactly where the data table contains a row for which the
            %value of the time column equals one of the Values.

            arguments
                obj (1, 1) MBExperimentTimeConfigurable
                data (1, 1) table
            end

            rowsToKeep = data.time == obj.Values(1);

            for valIdx = 2:length(obj.Values)
                % If there are multiple values, then each row that matches
                % any of the values should be kept. This is equivalent to
                % logical disjunction ("or").

                curRowsToKeep = data.time == obj.Values(valIdx);
                rowsToKeep = or(rowsToKeep, curRowsToKeep);
            end

            % Return the identified subset of the data rows.

            rows = rowsToKeep;
        end % findSubsetOfData

        function filenameString = get.FilenameString(obj)
            filenameString = "Time";
            for valIdx = 1:length(obj.Values)
                filenameString = ...
                    filenameString + "_" + num2str(obj.Values(valIdx));
            end            
        end

        function value = getSingleValue(obj)
            arguments
                obj (1, 1) MBExperimentTimeConfigurable ...
                    {mustHaveSingleValue}
            end

            value = obj.Values(1);
        end

        function configurables = multiply(obj)
            configurables = createArray(1, length(obj.Values), ...
                "MBExperimentTimeConfigurable");
            for valuesIdx = 1:length(obj.Values)             
                configurables(valuesIdx).Values = obj.Values(valuesIdx);               
            end
        end

        function possibilities = numberOfValues(obj)
            possibilities = length(obj.Values);
        end

        function sum = plus(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentTimeConfigurable
                obj2 (1, 1) MBExperimentTimeConfigurable
            end

            values1 = obj1.Values;
            values2 = obj2.Values;
            sum = MBExperimentTimeConfigurable;
            sum.Values = [values1 values2];
        end

        function obj = set.Values(obj, val)
            % For simplicity and clarity, we will remove duplicate times
            % and sort them in ascending order, as is the standard format
            % for SSIT time spans.
            
            obj.Values = unique(sort(val));            
        end
    end % Public instance methods

    methods (Static)
        function varName = getVarName()
            varName = "Time";
        end

        function varType = getVarType()
            varType = "double";
        end
    end
end

function mustHaveSingleValue(configurable)
arguments
    configurable (1, 1) MBExperimentTimeConfigurable
end

if configurable.numberOfValues > 1
    error("A time configurable " + configurable.FilenameString + ...
        " contained multiple values where exactly one is required.");
end

end