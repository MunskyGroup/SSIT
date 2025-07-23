classdef RandomSEDStrategy < AbstractSequentialExperimentDesignStrategy
    methods
        function round = apportionObservations(obj, round)
            arguments
                obj
                round (1, 1) SequentialExperimentRound
            end

            % Determine the maximum number of quanta that can be
            % apportioned to the experiments in this round. Then calculate,
            % for each quantum up to that total, a random integer no
            % greater than the number of experiments (i.e., assign each
            % quantum to an experiment). For each experiment, count the
            % number of corresponding assignments and multiply by the 
            % quantum to determine the total number of observations to be
            % assigned to that experiment's configuration.

            K = length(round.Experiments);
            N = idivide(int32(obj.ObservationsPerExperiment), ...
                int32(obj.ObservationQuantum));
            assignments = randi(K, 1, N);            
            for experimentIdx = 1:K
                curObservations = sum(assignments == experimentIdx) * ...
                    obj.ObservationQuantum;
                round.Experiments(experimentIdx).Configuration.NumberOfObservations = ...
                    curObservations;
            end             
        end % apportionObservations
    end
    methods (Access = protected)
        function round = designRoundInternal(obj, round)
            arguments
                obj
                round (1, 1) SequentialExperimentRound
            end
            disp(obj.ObservationsPerExperiment / obj.ObservationQuantum)
        end
    end
end