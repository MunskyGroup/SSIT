% Gene Expression System Experiment Class Definition
classdef GeneExprSysExperiment
    properties
        NCells (1,1) integer {mustBePositive}       
        Input (1,1) string {mustBeNonempty}
        NextTimeDelta (1,1) integer {mustBePositive}       
    end
end