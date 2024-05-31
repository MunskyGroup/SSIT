% D-optimal design minimizes the volume of the joint CI of all parameters,
% i.e. the determinant of the inverse of the FIM.
classdef FIMOptimalityCriterionD < FIMOptimalityCriterion
    methods
        function metric_func = getMetricFunction(~)
            metric_func = @(A)-det(inv(A));
        end
    end
end