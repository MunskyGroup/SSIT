classdef FIMOptimizedSEDStrategy < AbstractSequentialExperimentDesignStrategy
    methods (Static)
        function result = usesData(~)
            % FIM-optimized strategies need data.
            
            result = true;
        end
    end % Public static methods
    
    methods (Access = protected)
        function round = designRoundInternal(obj, round)
            arguments
                obj
                round (1, 1) SequentialExperimentRound
            end

            % Step 2 (of ancestral designRound).
            % Perform any required analysis specific to the 
            % design strategy.

            % Step 2a. Using the FSP approach, maximize the likelihood of
            % the data so far obtained.

            % Step 2b. Use Metropolis-Hastings sampling to estimate the
            % posterior distribution of the data.

            % Step 2c. Calculate a set of Fisher information matrices, one
            % for each set of M-H results, using parameter samples selected
            % uniformly from the M-H chain and for all experiment 
            % configurations.

            
        end % designRoundInternal
    end
end