classdef MetropolisHastingsAlgorithmResults
    %MetropolisHastingsAlgorithmResults encapsulates the results of running
    %the Metropolis-Hastings algorithm for the purpose of exploring
    %parameter space.

    properties
        Acceptance
        ParametersLogSpace
        Samples
        Value
    end
end