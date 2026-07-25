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

    methods
        function model = applyToModel(obj, model)
            %applyToModel takes an SSIT model and configures it according
            %to this configurable. SSIT allows a time span to have multiple
            %values, and the setter of the Values property already ensures
            %that all times are unique and in ascending order.
            
            arguments
                obj                 
                model (1, 1) SSIT
            end

            % We apply the time configurable to the model by setting the
            % tSpan property of the latter accordingly.

            model.tSpan = obj.Values;
        end % applyToModel

        function disp(obj)           
            disp(['Time = ' num2str(obj.Values)])            
        end

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