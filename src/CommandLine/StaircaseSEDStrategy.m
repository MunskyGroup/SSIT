classdef StaircaseSEDStrategy < AbstractSequentialExperimentDesignStrategy
    methods
        function round = apportionObservations(obj, round)
            arguments
                obj
                round (1, 1) SequentialExperimentRound
            end

            % We will allocate quanta to the experiments in a "staircase"
            % fashion, i.e. in the sequence 0, 1, 2, 3, ... quanta until
            % the end of the configurations or the maximum number of
            % observations is reached. Note that the sum of 1, ..., N is
            % N(N+1)/2 and must be less than ObservationsPerExperiment, so
            % we can solve for N easily.

            quanta = idivide(int32(obj.ObservationsPerExperiment), ...
                int32(obj.ObservationQuantum));
            N = floor((-1 + sqrt(1 + (8 * double(quanta)))) / 2);

            K = length(round.Experiments);

            % Experiments 1, ..., N+1 will contain 0, ..., N observations;
            % any remaining experiments will contain 0.

            for experimentIdx = 1:(N+1)
                round.Experiments(experimentIdx).Configuration.NumberOfObservations = ...
                    (experimentIdx - 1) * obj.ObservationQuantum;
            end
            for experimentIdx = (N+2):K
                round.Experiments(experimentIdx).Configuration.NumberOfObservations = ...
                    0;
            end
        end % apportionObservations
    end % Public methods
    
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