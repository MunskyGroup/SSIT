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

    methods (Abstract)
        model = applyToModel(obj, model)        
        varName = getVarName(obj)
        value = getSingleValue(obj)
        configurables = multiply(obj)
        possibilities = numberOfValues(obj)
    end

    methods (Abstract, Static)        
        varType = getVarType()
    end
end