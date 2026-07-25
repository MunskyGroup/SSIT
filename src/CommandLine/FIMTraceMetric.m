classdef FIMTraceMetric < AbstractFIMMetric
    % This metric maximizes the trace of the FIM, so to be used for
    % minimization routines, it needs to calculate the negative of the
    % trace.

    methods (Access = protected)
        function metricResult = calculateForFIMInternal(~, FIM)
            metricResult = -trace(FIM);
        end
    end % Protected methods
end