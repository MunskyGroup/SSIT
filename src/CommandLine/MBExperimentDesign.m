classdef MBExperimentDesign
    %MBExperimentDesign encapsulates how many observations an experimenter
    %should seek to obtain within a single experiment (or single round of a
    %sequence of experiments). 

    properties (Dependent)
        Configurations
        NumberOfConfigurations
        NumberOfObservations
    end

    properties
        % uint64 is used here to automatically require that the number of
        % observations is always a non-negative integer.
        
        ConfigurationToNumObservationsMap = ...
            configureDictionary("MBExperimentConfiguration", "uint64")
    end

    methods (Access = private)
        function [Obs, nonTimeConfigurations, times] = ...
                buildObservationMatrix(obj, nonTimeConfigurations, times)
            % We need to extract the non-time configurations and
            % measurement times from the map. We know that there will be no
            % more non-time configurations nor measurement times than the 
            % number of keys, though possibly fewer of either or both. This
            % is because each configuration can only have a single value 
            % for each of its configurables (including time), so there will
            % be at most as many unique values of each configurable as the 
            % number of configurations.
            %
            % Non-time configurations and/or times of interest can be
            % provided, none of which need to be included in any of the
            % configurations. If provided, they will be used in creating
            % that axis of the matrix, rather than their counterparts in
            % the observations map.

            arguments
                obj (1, 1) MBExperimentDesign
                nonTimeConfigurations (1, :) ...
                    MBExperimentNonTimeConfiguration = createArray(...
                        0, 0, "MBExperimentNonTimeConfiguration");
                times (1, :) double = []
            end

            k = obj.ConfigurationToNumObservationsMap.keys();

            % If we need to fill either the configurations or the times, we
            % need to enumerate over all unique values. For the times, we
            % can simply populate an array and then extract the unique
            % values; for the non-time configurations, the easiest way is
            % to make them keys of a "dummy" dictionary and then extract
            % the keys.

            fillConfigurations = isempty(nonTimeConfigurations);
            if fillConfigurations
                nonTimeConfigMap = configureDictionary(...
                    "MBExperimentNonTimeConfiguration", "double");
            end

            fillTimes = isempty(times);
            if fillTimes
                times = zeros(1, length(k));
            end

            if fillConfigurations || fillTimes
                for keyIdx = 1:length(k)
                    curConfig = k(keyIdx);

                    if fillConfigurations
                        nonTimeConfigMap(...
                            curConfig.NonTimeConfiguration) = 0;                        
                    end

                    if fillTimes
                        times(keyIdx) = ...
                            curConfig.TimeConfigurable.getSingleValue();
                    end
                end

                if fillConfigurations
                    nonTimeConfigurations = nonTimeConfigMap.keys();
                end
            end % [Fill something]
            
            times = unique(sort(times));

            Obs = zeros(length(nonTimeConfigurations), length(times));
        end % buildObservationMatrix
    end % Private methods

    methods
        function obj = excludeConfiguration(obj, configuration)
            % excludeConfiguration removes a configuration completely from
            % the design. This is NOT equivalent to setting the number of
            % corresponding observations to zero; the latter is a weaker
            % change, since the configuration will still be included in the
            % list of candidate configurations.

            arguments
                obj (1, 1) MBExperimentDesign
                configuration (1, 1) MBExperimentConfiguration ...
                    {mustHaveSingleValueForAllConfigurables(configuration)}
            end

            % If the configuration is not currently present in the map,
            % this will have no effect but will also raise no error.

            obj.ConfigurationToNumObservationsMap.remove(configuration);
        end % excludeConfiguration

        function configs = get.Configurations(obj)
            % MATLAB behaves strangely in returning a 0x1 array rather than
            % a 0x0 array when there are no keys. This causes issues in
            % downstream methods like MBExperimentConfiguration.disp:

            if obj.NumberOfConfigurations > 0
                configs = obj.ConfigurationToNumObservationsMap.keys();
            else
                configs = createArray(0, 0, "MBExperimentConfiguration");
            end
        end

        function numConfigs = get.NumberOfConfigurations(obj)
            numConfigs = obj.ConfigurationToNumObservationsMap.numEntries;
        end

        function numObservations = get.NumberOfObservations(obj)
            numObservations = ...
                sum(obj.ConfigurationToNumObservationsMap.values);
        end

        function [Obs, nonTimeConfigurations, times] = ...
                getAsObservationMatrix(obj, nonTimeConfigurations, times)
            % getAsObservationMatrix returns a matrix containing the number
            % of observations to be made at each non-time configuration 
            % (row) and measurement time (column). This is the format 
            % generated by methods such as optimizeObservationCounts.            
            %
            % The user can provide non-time configurations and/or times of
            % interest, none of which need to be included in any of the
            % configurations. If either is provided, it will be used in 
            % creating the relevant axis of the matrix, rather than its 
            % counterpart in the observations map.

            arguments
                obj (1, 1) MBExperimentDesign
                nonTimeConfigurations (1, :) ...
                    MBExperimentNonTimeConfiguration = createArray(...
                        0, 0, "MBExperimentNonTimeConfiguration");
                times (1, :) double = []
            end

            [Obs, nonTimeConfigurations, times] = ...
                obj.buildObservationMatrix(nonTimeConfigurations, times);
            Obs = uint64(Obs);

            for nonTimeConfigIdx = 1:length(nonTimeConfigurations)
                for timeIdx = 1:length(times)
                    curConfig = MBExperimentConfiguration();
                    curConfig.NonTimeConfiguration = ...
                        nonTimeConfigurations(nonTimeConfigIdx);
                    curConfig.TimeConfigurable.Values = times(timeIdx);

                    Obs(nonTimeConfigIdx, timeIdx) = ...
                        obj.ConfigurationToNumObservationsMap(curConfig);
                end
            end
        end % getAsObservationMatrix

        function numObservationsMap = ...
                getMostObservationsAtAnyTimeForNonTimeConfigurations(obj)
            arguments
                obj (1, 1) MBExperimentDesign                
            end

            numObservationsMap = configureDictionary(...
                "MBExperimentNonTimeConfiguration", "uint64");

            k = obj.ConfigurationToNumObservationsMap.keys();
            for keyIdx = 1:length(k)
                curKey = k(keyIdx);
                curNonTimeConfig = curKey.NonTimeConfiguration;
                if numObservationsMap.isKey(curNonTimeConfig)
                    curValue = numObservationsMap(curNonTimeConfig);
                    numObservationsMap(curNonTimeConfig) = max(...
                        curValue, ...
                        obj.ConfigurationToNumObservationsMap(curKey));
                else
                    numObservationsMap(curNonTimeConfig) = ...
                        obj.ConfigurationToNumObservationsMap(curKey);
                end
            end            
        end % getMostObservationsAtAnyTimeForNonTimeConfigurations

        function observations = getObservationsForConfiguration(...
                obj, configuration)
            arguments
                obj (1, 1) MBExperimentDesign
                configuration (1, 1) MBExperimentConfiguration ...
                    {mustHaveSingleValueForAllConfigurables(configuration)}                
            end

            if obj.ConfigurationToNumObservationsMap.isKey(configuration)
                observations = ...
                    obj.ConfigurationToNumObservationsMap(configuration);
            else
                observations = 0;
            end
        end % getObservationsForConfiguration

        function [Obs, nonTimeConfigurations, times] = ...
                getTheoreticalObservationMatrix(...
                obj, nonTimeConfigurations, times)
            % getTheoreticalObservationMatrix returns a matrix containing
            % the theoretical maximum number of observations that could be
            % made at each non-time configuration (row) and measurement 
            % time (column); this is the format generated by methods such 
            % as optimizeObservationCounts. This method exists because not 
            % all combinations of non-time configuration and measurement 
            % time might be available, so zero observations are 
            % theoretically possible in such cases, whereas infinitely many
            % observations are possible generally.
            %
            % The user can provide non-time configurations and/or times of
            % interest, none of which need to be included in any of the
            % configurations. If either is provided, it will be used in 
            % creating the relevant axis of the matrix, rather than its 
            % counterpart in the observations map.

            arguments
                obj (1, 1) MBExperimentDesign
                nonTimeConfigurations (1, :) ...
                    MBExperimentNonTimeConfiguration = createArray(...
                        0, 0, "MBExperimentNonTimeConfiguration");
                times (1, :) double = []
            end

            % buildObservationMatrix returns an m-by-n matrix of zeros.

            [Obs, nonTimeConfigurations, times] = ...
                obj.buildObservationMatrix(nonTimeConfigurations, times);
            
            for nonTimeConfigIdx = 1:length(nonTimeConfigurations)
                for timeIdx = 1:length(times)
                    curConfig = MBExperimentConfiguration();
                    curConfig.NonTimeConfiguration = ...
                        nonTimeConfigurations(nonTimeConfigIdx);
                    curConfig.TimeConfigurable.Values = times(timeIdx);

                    % If the current configuration IS a key of the map,
                    % that means it is an available configuration, and
                    % therefore, arbitrarily many observations could be
                    % possible. If it is NOT a key, then no observations
                    % are possible, so we will leave the corresponding
                    % entry of the observation matrix at zero.

                    if obj.ConfigurationToNumObservationsMap.isKey(curConfig)
                        Obs(nonTimeConfigIdx, timeIdx) = inf;
                    end
                end
            end
        end % getTheoreticalObservationMatrix

        function obj = MBExperimentDesign(configurations)
            arguments
                configurations (1, :) MBExperimentConfiguration = ...
                    createArray(0, 0, "MBExperimentConfiguration")
            end

            % Create an "empty" design, with zero observations for every
            % supplied configuration (if any).

            for configIdx = 1:length(configurations)
                curConfig = configurations(configIdx);
                obj.ConfigurationToNumObservationsMap(curConfig) = 0;
            end
        end % MBExperimentDesign

        function sum = plus(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentDesign
                obj2 (1, 1) MBExperimentDesign
            end

            % The addition of two experiment designs is equivalent to
            % merging their dictionaries, i.e., taking the union of all
            % keys and using values that are the sums of the values in the
            % dictionaries being merged.

            map1 = obj1.ConfigurationToNumObservationsMap;
            map2 = obj2.ConfigurationToNumObservationsMap;

            sumMap = dictionary;

            keys1 = map1.keys();

            for keyIdx = 1:length(keys1)
                curKey = keys1(keyIdx);
                sumMap(curKey) = map1(curKey);
                if map2.isKey(curKey)
                    sumMap(curKey) = sumMap(curKey) + map2(curKey);
                end
            end

            % At this point, sumMap includes the correct value for all keys
            % that are in map1 (whether or not they are also in map2). We
            % therefore need to add the value for all keys that are only in
            % map2:

            keys2 = map2.keys();

            for keyIdx = 1:length(keys2)
                curKey = keys2(keyIdx);
                if ~map1.isKey(curKey) % Keys in map1 were already included 
                    sumMap(curKey) = map2(curKey);
                end                
            end

            % Now the maps have been merged.

            sum = MBExperimentDesign();
            sum.ConfigurationToNumObservationsMap = sumMap;
        end % plus

        function obj = setFromObservationMatrix(...
                obj, Obs, nonTimeConfigurations, times)
            %setFromObservationMatrix accepts a matrix containing the
            %number of observations to be made at each non-time 
            %configuration (row) and measurement time (column); this is the
            %format generated by methods such as optimizeObservationCounts.
            %It also requires a list of lists of non-time configurations
            %and measurement times corresponding to the matrix.

            arguments
                obj (1, 1) MBExperimentDesign
                Obs (:, :) uint64
                nonTimeConfigurations (1, :) ...
                    MBExperimentNonTimeConfiguration {mustBeNonempty} 
                times (1, :) double {mustBeNonempty, mustBeNonnegative}
            end

            [rows, cols] = size(Obs);
            if rows ~= length(nonTimeConfigurations)
                error("Cannot set design from matrix with " + ...
                    "incompatible number of non-time configurations")
            end
            if cols ~= length(times)
                error("Cannot set design from matrix with " + ...
                    "incompatible number of measurement times")
            end

            % Having verified that the inputs are mutually consistent,
            % reset the dictionary and then iterate over the observation
            % matrix, populating the dictionary with the corresponding
            % (key, value) pair.

            obj.ConfigurationToNumObservationsMap = ...
                configureDictionary("MBExperimentConfiguration", "uint64");

            for nonTimeConfigIdx = 1:rows
                for timeIdx = 1:cols
                    curConfig = MBExperimentConfiguration();
                    curConfig.NonTimeConfiguration = ...
                        nonTimeConfigurations(nonTimeConfigIdx);
                    curConfig.TimeConfigurable.Values = times(timeIdx);

                    obj.ConfigurationToNumObservationsMap(curConfig) = ...
                        Obs(nonTimeConfigIdx, timeIdx);
                end
            end
        end % setFromObservationMatrix

        function obj = setObservationsForConfiguration(...
                obj, configuration, observations)
            arguments
                obj (1, 1) MBExperimentDesign
                configuration (1, 1) MBExperimentConfiguration ...
                    {mustHaveSingleValueForAllConfigurables(configuration)}
                observations (1, 1) uint64
            end

            obj.ConfigurationToNumObservationsMap(configuration) = ...
                observations;
        end % setObservationsForConfiguration
    end % Public methods
end % MBExperimentDesign

function mustHaveSingleValueForAllConfigurables(configurations)
    arguments
        configurations (1, :) MBExperimentConfiguration
    end
    
    % Here, we validate that every configurable in each configuration only has 
    % a single value, via the getSingleValue method on each configurable, since
    % that method will error and print the values if there are multiple of
    % them.
    
    for configIdx = 1:length(configurations)
        curConfig = configurations(configIdx);
        nonTimeConfig = curConfig.NonTimeConfiguration;
        for nonTimeConfigIdx = 1:length(nonTimeConfig.NonTimeConfigurables)
            curNonTimeConfigurable = ...
                nonTimeConfig.NonTimeConfigurables(nonTimeConfigIdx);
            curNonTimeConfigurable.getSingleValue();
        end
    
        curConfig.TimeConfigurable.getSingleValue();
    end
end