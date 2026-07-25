classdef MetropolisHastingsAlgorithmOptions
    %MetropolisHastingsAlgorithmOptions encapsulates some options needed by
    %the algorithm runner.

    properties
        Progress (1, 1) logical = false;
        SaveFile (1, 1) string = "TMPMHChain_1.mat"
    end
end