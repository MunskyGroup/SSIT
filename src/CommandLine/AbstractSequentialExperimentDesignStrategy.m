classdef AbstractSequentialExperimentDesignStrategy
    properties
        ObservationsPerExperiment (1, 1) double {mustBeInteger, mustBePositive} = 100
        ObservationQuantum (1, 1) double {mustBeInteger, mustBePositive} = 10
    end
    methods
        function round = designRound(obj, round)            
            % Step 1. Perform any required analysis specific to the 
            % design strategy, including (if appropriate) the importing and
            % fitting of all relevant data collected thus far, i.e. all
            % data pertaining to the configuration of each experiment in 
            % the round.
            
            round = obj.designRoundInternal(round);

            % Step 2. Apportion the available/desired number of
            % observations across the experiments in the round.

            round = obj.apportionObservations(round);
        end % designRound        
    end % Public object methods

    methods (Static)
        function result = usesData(~)
            % Design strategies do not use data by default; descendants may
            % override this.

            result = false;
        end
    end % Public static methods

    methods (Abstract)
        round = apportionObservations(obj, round)        
    end
    
    methods (Access = protected)
        function round = designRoundInternal(obj, round)
            arguments
                obj
                round (1, 1) SequentialExperimentRound
            end

            % In this ultimate ancestor, the only thing we will do is to
            % load data into the models in the round IF the strategy uses
            % data.

            if obj.usesData
            end
        end
    end
end