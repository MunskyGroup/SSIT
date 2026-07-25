classdef FIMEMetric < AbstractFIMMetric
    % This metric maximizes the smallest eigenvalue of the FIM, so to be
    % used for minimization routines, it needs to calculate the negative of
    % the minimum of the eigenvalues.

    methods (Access = protected)
        function metricResult = calculateForFIMInternal(~, FIM)
            metricResult = -min(eig(FIM));
        end
    end % Protected methods
end