classdef FIMDMetric < AbstractFIMMetric
    % This metric maximizes the expected determinant of the FIM, so to be
    % used for minimization routines, it needs to calculate the negative of
    % the determinant.

    methods (Access = protected)
        function metricResult = calculateForFIMInternal(~, FIM)
            metricResult = -max(0, det(FIM));
        end
    end % Protected methods
end