classdef AbstractExperimentConfigurable < matlab.mixin.Heterogeneous
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