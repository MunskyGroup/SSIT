% A-optimal design minimizes total parameter variance, i.e. the trace of
% the inverse of the FIM.
classdef FIMOptimalityCriterionA < FIMOptimalityCriterion
    methods
        function metric_func = getMetricFunction(~)
            metric_func = @(A)-trace(inv(A));
        end
    end
end