classdef MBExperimentConfiguration
    %MBExperimentConfiguration encapsulates a configuration, i.e. an
    %assignment of a single value for all "configurables": the non-time
    %configurations (input expressions and designed parameters of a model)
    %and the measurement time(s).
    
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
                data (1, 1) table
            end

            rowsToKeep = obj.TimeConfigurable.findSubsetOfData(data);

            nonTimeConfigs = obj.NonTimeConfiguration.NonTimeConfigurables;
            for configIdx = 1:length(nonTimeConfigs)
                % All rows must match, i.e., have the desired value for
                % each non-time configurable and one of the desired values
                % for the time configurable. This is equivalent to logical
                % conjunction ("and").
                
                curRowsToKeep = ...
                    nonTimeConfigs(configIdx).findSubsetOfData(data);
                rowsToKeep = and(rowsToKeep, curRowsToKeep);
            end

            % Return the identified subset of the data rows and all data
            % columns.

            data = data(rowsToKeep, :);
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
    end
end