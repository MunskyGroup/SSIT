classdef MBExperimentConfiguration
    %MBExperimentConfiguration encapsulates a configuration, i.e. an
    %assignment of a single value for all "configurables": the non-time
    %configurations (input expressions and designed parameters of a model)
    %and the measurement time(s).

    properties (Dependent)
        FilenameString
    end
    
    properties
        NonTimeConfiguration (1, 1) MBExperimentNonTimeConfiguration
        TimeConfigurable (1, 1) MBExperimentTimeConfigurable
    end

    methods
        function data = applyToData(obj, data)
            %applyToData takes a dataset in the form of a MATLAB table and
            %subsets it according to this configuration; in other words,
            %the output table will only contain those rows for which
            %   1. The value of each column corresponding to a non-time
            %   configurable will equal the single value of that
            %   configurable.
            %   AND
            %   2. The value of each column corresponding to the time
            %   configurable will equal one of the values of the time
            %   configurable.

            arguments
                obj (1, 1) MBExperimentConfiguration
                data table
            end

            % Identify and return the appropriate subset of the data rows 
            % and all data columns.

            if ~isempty(data)
                data = data(obj.findSubsetOfData(data), :);
            end
        end % applyToData

        function model = applyToModel(obj, model)
            %applyToModel takes a model and configures it according to
            %these configurables.
            arguments
                obj (1, 1) MBExperimentConfiguration              
                model (1, 1) MBExperimentModel
            end

            nonTimeConfigurables = ...
                obj.NonTimeConfiguration.NonTimeConfigurables;
            for nonTimeConfigIdx = 1:length(nonTimeConfigurables)
                curConfig = nonTimeConfigurables(nonTimeConfigIdx);
                model = curConfig.applyToModel(model);
            end

            model = obj.TimeConfigurable.applyToModel(model);
        end % applyToModel

        function disp(obj)
            disp(obj.NonTimeConfiguration)
            disp(obj.TimeConfigurable)
        end % disp

        function isEqual = eq(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentConfiguration
                obj2 (1, 1) MBExperimentConfiguration
            end

            isEqual = ...
                (obj1.NonTimeConfiguration == ...
                obj2.NonTimeConfiguration) && ...
                (obj1.TimeConfigurable == obj2.TimeConfigurable);
        end

        function rows = findSubsetOfData(obj, data)
            %findSubsetOfData takes a dataset in the form of a MATLAB table
            %and finds its subset that matches this configuration:
            %
            %   1. The value of each column corresponding to a non-time
            %   configurable will equal the single value of that
            %   configurable.
            %   AND
            %   2. The value of each column corresponding to the time
            %   configurable will equal one of the values of the time
            %   configurable.

            arguments
                obj (1, 1) MBExperimentConfiguration
                data table
            end

            rows = [];
            if ~isempty(data)
                rows = obj.TimeConfigurable.findSubsetOfData(data);
    
                nonTimeConfigs = ...
                    obj.NonTimeConfiguration.NonTimeConfigurables;
                for configIdx = 1:length(nonTimeConfigs)
                    % All rows must match, i.e., have the desired value for
                    % each non-time configurable and one of the desired values
                    % for the time configurable. This is equivalent to logical
                    % conjunction ("and").
    
                    curRowsToKeep = ...
                        nonTimeConfigs(configIdx).findSubsetOfData(data);
                    rows = and(rows, curRowsToKeep);
                end
            end % [Non-empty data]
        end % findSubsetOfData

        function filenameString = get.FilenameString(obj)
            %FilenameString represents the configuration as a string
            %suitable for use in a filename.
           
            filenameString = join([...
                obj.NonTimeConfiguration.FilenameString ...
                obj.TimeConfigurable.FilenameString], ...
                "_");
        end

        function obj = setConfigurable(obj, configurable)
            arguments
                obj (1, 1) MBExperimentConfiguration
                configurable (1, 1) MBAbstractExperimentConfigurable
            end

            if isa(configurable, "MBExperimentTimeConfigurable")
                obj.TimeConfigurable = configurable;
            else
                obj.NonTimeConfiguration.NonTimeConfigurables(end + 1) ...
                    = configurable;
            end
        end % setConfigurable
    end % Public methods
end