classdef RandomSEDStrategy < AbstractSequentialExperimentDesignStrategy
    methods
        function round = apportionObservations(obj, round)
            arguments
                obj
                round (1, 1) SequentialExperimentRound
            end

            % Determine the maximum number of quanta that can be
            % apportioned to each experiment. Then calculate, for each
            % experiment, a random integer no greater than that maximum.

            K = length(round.Experiments);
            N = idivide(int32(obj.ObservationsPerExperiment), ...
                int32(obj.ObservationQuantum));
            quanta = randi(N, 1, K);
            observations = quanta * obj.ObservationQuantum;
            for experimentIdx = 1:K
                round.Experiments(experimentIdx).Configuration.NumberOfObservations = ...
                    observations(experimentIdx);
            end             
        end
    end    
end