classdef AbstractSequentialExperimentDesignStrategy
    properties
        ObservationsPerExperiment (1, 1) double {mustBeInteger, mustBePositive} = 100
        ObservationQuantum (1, 1) double {mustBeInteger, mustBePositive} = 10
    end
    methods
        function round = designRound(obj, round)
            % Step 1. Import all relevant data collected thus far, i.e.
            % all data pertaining to the configuration of each experiment
            % in the round.

            % Step 2. Perform any required analysis specific to the 
            % design strategy.
            
            round = obj.designRoundInternal(round);

            % Step 3. Apportion the available/desired number of
            % observations across the experiments in the round.

            round = obj.apportionObservations(round);
        end
    end
    methods (Abstract)
        round = apportionObservations(obj, round)        
    end
    methods (Access = protected)
        function round = designRoundInternal(~, round)
            arguments
                ~
                round (1, 1) SequentialExperimentRound
            end
        end
    end
end