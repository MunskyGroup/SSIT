classdef AbstractExperimentConfigurable < matlab.mixin.Heterogeneous
    methods (Abstract)
        configurables = multiply(obj)
        possibilities = numberOfValues(obj)
    end
end