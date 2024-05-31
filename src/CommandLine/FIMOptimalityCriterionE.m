% E-optimal design minimizes the largest eigenvalue of the FIM, thereby
% minimizing the uncertainties in the worst-case direction in parameter
% space.
classdef FIMOptimalityCriterionE < FIMOptimalityCriterion
    methods
        function metric_func = getMetricFunction(~)
            metric_func = @(A)-max(eig(A));
        end
    end
end