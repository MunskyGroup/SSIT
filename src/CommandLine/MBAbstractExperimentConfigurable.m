classdef MBAbstractExperimentConfigurable < matlab.mixin.Heterogeneous
    %MBAbstractExperimentConfigurable is an ancestor of the various types
    %of experiment configurables, i.e., those aspects of the system and
    %experiment that are within the control of an experimenter.
    %
    %   Two types of configurables are supported:
    %   - Non-time configurables (input expressions and parameters)
    %   - Measurement time(s)

    properties (Abstract, Dependent)
        FilenameString
    end

    methods (Abstract, Access = protected)
        isEqual = equalsSibling(obj1, obj2)
    end

    methods (Abstract)        
        model = applyToModel(obj, model)       
        rows = findSubsetOfData(obj, data)
        varName = getVarName(obj)
        value = getSingleValue(obj)
        configurables = multiply(obj)
        possibilities = numberOfValues(obj)
    end

    methods (Abstract, Static)        
        varType = getVarType()
    end

    methods
        function isEqual = eq(obj1, obj2)
            arguments
                obj1 (1, 1) MBAbstractExperimentConfigurable
                obj2 (1, 1) MBAbstractExperimentConfigurable
            end

            % First, verify that the two objects are of exactly the same
            % class; if so, then call that class' equalsSibling method to
            % determine the final answer. We achieve this by
            % short-circuiting the conjunction:

            isEqual = strcmp(class(obj1), class(obj2)) && ...
                equalsSibling(obj1, obj2);
        end % eq
    end % Public methods
end