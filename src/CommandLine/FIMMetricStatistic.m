classdef FIMMetricStatistic
    enumeration
        Mean
        Median
    end
    
    methods
        function stat = calculateStatisticOfObjective(...
                obj, objectiveFunction)
            arguments
                obj (1, 1) FIMMetricStatistic
                objectiveFunction {mustBeMatrix}
            end
            switch obj
                case "Mean"
                    stat = min(mean(objectiveFunction, 2));
                case "Median"
                    stat = min(median(objectiveFunction, 2));
                otherwise
                    error("Statistic calculation not supported for " + obj)
            end
        end
    end
end