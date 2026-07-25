classdef FIMSubsetInverseDMetric < FIMCovarianceDMetric & ...
        AbstractFIMSubsetMetric
    % This metric minimizes the expected determinant of the MLE covariance,
    % so to be used for minimization routines, it needs to calculate the
    % determinant of the inverse of the FIM.

    methods (Access = protected)
        function metricResult = calculateForFIMInternal(~, FIM)
            metricResult = max(0, det(inv(FIM)));
        end
    end % Protected methods
end