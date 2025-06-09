classdef SequentialExperimentDesigner
    properties
        DataType (1, 1) ExperimentalDataType = Simulated
        Rounds (1, :) SequentialExperimentRound
        Strategy (1, 1) AbstractSequentialExperimentDesignStrategy = RandomSEDStrategy()
    end

    methods
        function round = designNextRound(obj)
            round = SequentialExperimentRound()
            
        end

        function performNextRound(obj)
        end
    end
end