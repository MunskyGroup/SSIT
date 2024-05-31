classdef FIMOptimalityCriterionUnclear < FIMOptimalityCriterion
    properties (SetAccess = private)
        Desc (1,1) string {mustBeNonempty}
        NumTimepoints (1,1) int {mustBePositive}
    end
    
    methods
        % Need constructor to allow for extra properties
        function c = FIMOptimalityCriterionUnclear(desc, nT)
            % Call superclass constructor
            c@FIMOptimalityCriterion({});

            c.Desc = desc;
            c.NumTimepoints = nT;
        end

        function metric_func = getMetricFunction(criterion)
            if strcmp(criterion.Desc(1:2),'TR')
                k = eval(criterion.Desc(3:end));
                metric_func = @(A)det(inv(A(k,k)));
            else  % all parameters are free.
                k = eval(criterion.Desc);
                ek = zeros(length(k), criterion.NumTimepoints);
                ek(1:length(k),k) = eye(length(k));
                metric_func = @(A)det(ek*inv(A)*ek');
            end
        end
    end
end