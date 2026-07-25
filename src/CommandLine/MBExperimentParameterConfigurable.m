classdef MBExperimentParameterConfigurable < ...
        MBAbstractExperimentNonTimeConfigurable
    %MBExperimentParameterConfigurable represents an individual parameter
    %configurable, i.e., a kinetic parameter with a specific name and one
    %or more possible values.
    %
    %   Example:
    %       InputName = "t_tpl"
    %       Values = [-10 10 30.5 60 75]
    %   In this example, the "t_tpl" model parameter can take one of five
    %   possible values. Observe that the values need not be integral nor
    %   evenly spaced; however, they must be numeric.

    properties
        ParameterName (1, 1) string {mustBeNonempty} = "ParameterName"
        Values (1, :) double {mustBeNonempty} = [0]
    end

    properties (Dependent)
        FilenameString
    end

    methods
        function model = applyToModel(obj, model)
            %applyToModel takes an SSIT model and configures it according
            %to this configurable, provided that it has only a single
            %value, since SSIT does not allow a single parameter to take
            %multiple values.           
            
            arguments
                obj (1, 1) MBExperimentParameterConfigurable ...
                    {mustHaveSingleValue}                
                model (1, 1) SSIT
            end

            value = obj.getSingleValue();

            % We apply the parameter to the model by setting the
            % corresponding value in its "parameters" cell array.

            parameterNameMatchIndices = strcmp(model.parameters(:, 1), ...
                obj.ParameterName);
            matchIdx = find(parameterNameMatchIndices, 1);
            if ~isempty(matchIdx)
                model.parameters(:, 2) = value;
            else
                % The parameter does not yet exist in the model. Append it
                % and its value to the cell array.
                [numberOfParameters, ~] = size(pars);
                model.parameters{numberOfParameters + 1, 1} = ...
                    obj.ParameterName;
                model.parameters{numberOfParameters + 1, 2} = ...
                    value;
                matchIdx = numberOfParameters + 1;
            end

            % We also ensure that the parameter is not included in the list
            % of parameters to be fit.

            parametersToFit = model.fittingOptions.modelVarsToFit;
            if strcmp(parametersToFit, 'all') || isempty(parametersToFit)
                % Obtain a list of all parameters in the model.
                parametersToFit = 1:size(model.parameters, 1);
            end

            % Transpose the parameter vector to a row vector and remove the
            % row (if it exists) containing the configured parameter. Then
            % update the modelVarsToFit property with the trimmed vector in
            % column form.

            parametersToFit = parametersToFit';
            parametersToFitIdx = find(parametersToFit == matchIdx);
            parametersToFit = ...
                removerows(parametersToFit, parametersToFitIdx);
            model.fittingOptions.modelVarsToFit = parametersToFit';
        end % applyToModel

        function disp(obj)
            fprintf('%s = %s\n', obj.ParameterName, obj.Values)           
        end

        function filenameString = get.FilenameString(obj)
            %FilenameString represents the configurable as a string
            %suitable for use in a filename.

            filenameString = join([obj.ParameterName obj.Values], "_");
        end

        function value = getSingleValue(obj)
            arguments
                obj (1, 1) MBExperimentParameterConfigurable ...
                    {mustHaveSingleValue}
            end
            
            value = obj.Values(1);        
        end

        function varName = getVarName(obj)
            varName = obj.ParameterName;
        end

        function result = lt(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentParameterConfigurable
                obj2 (1, 1) MBExperimentParameterConfigurable
            end

            if obj1.ParameterName < obj2.InputName
                result = true;
            elseif obj1.ParameterName > obj2.InputName
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
            %containing a distinct single value in the original parameter
            %configurable.
            configurables = createArray(1, length(obj.Values), ...
                "MBExperimentParameterConfigurable");
            for valuesIdx = 1:length(obj.Values)
                configurables(valuesIdx).InputName = obj.ParameterName;
                configurables(valuesIdx).Values = obj.Values(valuesIdx);               
            end
        end

        function possibilities = numberOfValues(obj)            
            possibilities = length(obj.Values);
        end

        function obj = set.Values(obj, val)
            % For simplicity and clarity, we will remove duplicate
            % parameter values and sort them in ascending order.

            obj.Values = unique(sort(val));            
        end
    end % Public instance methods

    methods (Static)
        function varType = getVarType()
            varType = "double";
        end
    end
end

function mustHaveSingleValue(configurable)
arguments
    configurable (1, 1) MBExperimentParameterConfigurable
end

if configurable.numberOfValues > 1
    error("A parameter configurable " + configurable.FilenameString + ...
        " contained multiple values where exactly one is required.");
end

end