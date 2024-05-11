% Gene Expression System Experiment Class Definition
classdef GeneExprSysExperiment
    properties
        NCells (1,1) integer {mustBePositive}       
        InputIdx (1,1) integer {mustBePositive} % MATLAB indices start at 1
        NextTimeDelta (1,1) integer {mustBePositive}       
    end
end