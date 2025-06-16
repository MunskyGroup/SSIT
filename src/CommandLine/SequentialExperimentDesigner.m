classdef SequentialExperimentDesigner
    properties        
        DataType (1, 1) ExperimentalDataType = ExperimentalDataType.Simulated
        ExperimentalConfigurations (1, :) ExperimentConfiguration
        IsTesting (1, 1) logical = false
        MaximumNumberOfRounds (1, 1) double {mustBeInteger, mustBePositive} = intmax("uint64")
        MaximumNumberOfObservations (1, 1) double {mustBeInteger, mustBePositive} = intmax("uint64")
        Model (1, 1) DiscoverableModel        
        RNGSeed (1, 1) double {mustBeInteger} = 42
        Strategy (1, 1) AbstractSequentialExperimentDesignStrategy = RandomSEDStrategy()
    end
    properties (SetAccess = private)
        DataLocations (1, :) string = []
        DataStore %(1, 1) = datastore()
        NumberOfObservationsCompleted (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
        NumberOfRoundsCompleted (1, 1) double {mustBeInteger, mustBeNonnegative} = 0
        Rounds (1, :) SequentialExperimentRound
    end

    methods
        function [obj, round] = designNextRound(obj)
            if obj.NumberOfRoundsCompleted == obj.MaximumNumberOfRounds
                disp('Maximum number of experiment rounds completed')
            else
                if obj.NumberOfRoundsCompleted == length(obj.Rounds)
                    % We need to grow the array of Rounds. To preserve all
                    % existing results, simply concatenate the existing
                    % Rounds array with a new template round.
                    obj.Rounds = [obj.Rounds ...
                        SequentialExperimentRound(obj.ExperimentalConfigurations)];
                end
                    round = obj.Strategy.designRound(...
                        obj.Rounds(1 + obj.NumberOfRoundsCompleted));
                    obj.Rounds(1 + obj.NumberOfRoundsCompleted) = round;
            end            
        end

        function obj = performNextRound(obj, location)
            arguments
                obj 
                location (1, :) string {mustBeNonempty}
            end      

            roundPerformed = false;

            switch obj.DataType
                case Empirical
                    % When working with empirical data, the caller must
                    % supply the location(s) of the data file(s) that will
                    % have been collected according to the design of this
                    % round. Here, performing the round simply means
                    % incorporating the new data into the datastore.
                    % Because we always create a new datastore, using the
                    % new set of all locations, the datastore is always an
                    % initial, fully readable state for the design of the
                    % next round, enabling immediate loading and fitting of
                    % all data.

                    if ~isempty(location)
                        obj.DataLocations = [obj.DataLocations location];
                        obj.DataStore = datastore(obj.DataLocations);
                        obj.NumberOfRoundsCompleted = 1 + ...
                            obj.NumberOfRoundsCompleted;
                        roundPerformed = true;
                    else
                        disp('performNextRound requires at least one data location when working with empirical data!')
                    end
                case Simulated
                    % TO DO: Run SSA simulations according to the true
                    % model, write these to files, and then append them to
                    % the datastore.
            end % switch obj.DataType

            if roundPerformed
                obj.NumberOfObservationsCompleted = ...
                    obj.NumberOfObservationsCompleted + ...
                    obj.Rounds(obj.NumberOfRoundsCompleted).NumberOfObservations;
            end
        end % performNextRound
    end % Public methods
end