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
    methods
        function likelihood = getLikelihood(obj, x, hasData, stateSpace)
            arguments
                obj
                x
                hasData (1, 1) logical
                stateSpace = []
            end

            likelihood = 0;
            includePrior = true;
            if hasData
                if includePrior
                    likelihood = likelihood + obj.computeLikelihood(...
                        exp(x), stateSpace);
                    includePrior = false;
                else
                    tmpModel = obj;
                    tmpModel.fittingOptions.logPrior = [];
                    likelihood = likelihood + tmpModel.computeLikelihood(...
                        exp(x), stateSpace);
                end
            end
        end % getLikelihood
    end % methods
end