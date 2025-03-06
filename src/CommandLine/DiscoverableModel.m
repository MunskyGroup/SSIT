classdef DiscoverableModel < SSIT
    properties
        ArchitectureName (1,1) string
        DataToFit (:,2) cell
        FitParameters (1,:) double {mustBePositive, mustBeInteger}
        MuLog10Prior (1,:) double
        % Number of Metropolis-Hastings Samples to run
        NumberOfMHSamples (1,1) double {mustBePositive, mustBeInteger} 
        NumberOfTimepoints (1,1) double {mustBePositive, mustBeInteger}
        SigmaLog10Prior (1,:) double {mustBeNonnegative}
        % TrueParameters (:, 2) cell (UNNEEDED?)
    end
end