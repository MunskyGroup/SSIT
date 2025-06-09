classdef SequentialExperimentRound
    properties
        Experiments (1, :) DesignedExperiment
    end
    properties %(Access = private)
        NumberOfConfigs (1, 1) double {mustBeInteger, mustBePositive} = 1
        NumberOfConfigurables (1, 1) double {mustBeInteger, mustBePositive} = 1
        VarNames (1, :) string {mustBeNonempty} = 'VarNamePlaceholder'
        VarTypes (1, :) string {mustBeNonempty} = 'VarTypePlaceholder'
    end
    methods
        function obj = SequentialExperimentRound(configs)
            arguments
                configs (1, :) ExperimentConfiguration
            end

            obj.NumberOfConfigs = length(configs);
            obj.NumberOfConfigurables = length(configs(1).Configurables);
            
            obj.VarTypes = ...
                createArray(1, obj.NumberOfConfigurables + 1, "string");
            obj.VarNames = ...
                createArray(1, obj.NumberOfConfigurables + 1, "string");
            obj.VarTypes(end) = 'int64';
            obj.VarNames(end) = 'Number of Observations';

            for configurableIdx = 1:obj.NumberOfConfigurables
                obj.VarTypes(configurableIdx) = ...
                    configs(1).Configurables(configurableIdx).getVarType;
                obj.VarNames(configurableIdx) = ...
                    configs(1).Configurables(configurableIdx).getVarName;
            end

            obj.Experiments = ...
                createArray(1, obj.NumberOfConfigs, "DesignedExperiment");
            for configIdx = 1:obj.NumberOfConfigs
                obj.Experiments(configIdx).Configuration = ...
                    configs(configIdx);
            end
        end
        function t = exportAsTable(obj)
            t = table(...
                'Size', [obj.NumberOfConfigs (obj.NumberOfConfigurables + 1)], ...
                'VariableTypes', obj.VarTypes);
            t.Properties.VariableNames = obj.VarNames;
            for configIdx = 1:obj.NumberOfConfigs
                for configurableIdx = 1:obj.NumberOfConfigurables
                    t(configIdx, configurableIdx) = ...
                        {obj.Experiments(configIdx).Configuration.Configurables(configurableIdx).getSingleValue};
                end
                t(configIdx, end) = {obj.Experiments(configIdx).Configuration.NumberOfObservations};
            end
        end
    end
end