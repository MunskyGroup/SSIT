classdef MBExperimentNonTimeConfiguration
    %MBExperimentNonTimeConfiguration encapsulates the non-time-specific
    %portions of an experiment configuration, i.e. an assignment of single
    %values for all designed parameters and for all input expressions of a 
    %model.

    properties
        NonTimeConfigurables (1, :) MBAbstractExperimentNonTimeConfigurable
    end
end