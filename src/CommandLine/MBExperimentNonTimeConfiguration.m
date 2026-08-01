classdef MBExperimentNonTimeConfiguration
    %MBExperimentNonTimeConfiguration encapsulates the non-time-specific
    %portions of an experiment configuration, i.e. an assignment of single
    %values for all designed parameters and for all input expressions of a 
    %model.

    properties (Dependent)
        FilenameString
    end

    properties
        NonTimeConfigurables (1, :) ...
            MBAbstractExperimentNonTimeConfigurable ...
            {mustBeValidConfiguration}
    end

    methods
        function disp(obj)
            if ~isempty(obj.NonTimeConfigurables)
                disp(obj.NonTimeConfigurables(1))
                for configIdx = 2:length(obj.NonTimeConfigurables)                    
                    disp(obj.NonTimeConfigurables(configIdx))
                end
            end            
        end % disp

        function isEqual = eq(obj1, obj2)
            arguments
                obj1 (1, 1) MBExperimentNonTimeConfiguration
                obj2 (1, 1) MBExperimentNonTimeConfiguration
            end

            isEqual = true;

            if length(obj1.NonTimeConfigurables) ~= ...
                    length(obj2.NonTimeConfigurables)
                isEqual = false;
            else
                for config1Idx = 1:length(obj1.NonTimeConfigurables)
                    curConfig1 = obj1.NonTimeConfigurables(config1Idx);
                    foundMatch = false;

                    for config2Idx = 1:length(obj2.NonTimeConfigurables)
                        curConfig2 = obj2.NonTimeConfigurables(config2Idx);
                        if curConfig1 == curConfig2
                            foundMatch = true;
                            break
                        end
                    end

                    if ~foundMatch
                        isEqual = false;
                        break
                    end
                end
            end % [Lengths match]
        end % eq (==)

        function filenameString = get.FilenameString(obj)
            %FilenameString represents the configuration as a string
            %suitable for use in a filename.

            configurableStrings = ...
                createArray(1, length(obj.NonTimeConfigurables), "string");
            for configIdx = 1:length(obj.NonTimeConfigurables)
                configurableStrings(configIdx) = ...
                    obj.NonTimeConfigurables(configIdx).FilenameString;
            end

            filenameString = join(configurableStrings, "_");
        end
    end % Public methods
end % class

function mustBeValidConfiguration(configurables)
arguments
    configurables (1, :) MBAbstractExperimentNonTimeConfigurable
end

% To be a valid configuration, a list of non-time configurations must
% satisfy the following properties:
%   1. Each configurable has a single value.
%   2. No two configurables have the same name. Even if different types of
%   configurables (e.g., an input expression and designed parameter) shared
%   the same name, collisions would result when trying to use the model.

numConfigs = length(configurables);
configNames = createArray(1, numConfigs, "string");

for configIdx = 1:numConfigs
    curConfig = configurables(configIdx);
    curConfig.getSingleValue(); % This will error if multiple values exist.
    configNames(configIdx) = curConfig.getVarName;
end

uniqueConfigNames = unique(configNames);
if length(uniqueConfigNames) < numConfigs
    error("A non-time configuration had duplicate configurable " + ...
        "names: " + join(configNames, ", "))
end

end