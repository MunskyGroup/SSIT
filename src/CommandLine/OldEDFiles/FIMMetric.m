classdef FIMMetric
    properties
        MetricKind (1,1) FIMMetricKind
        Indices (1,:) {mustBePositive, mustBeInteger} = []
    end
end