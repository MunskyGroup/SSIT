classdef FIMSubsetDMetric < FIMDMetric & AbstractFIMSubsetMetric
    % This metric maximizes the expected determinant of the FIM, so to be
    % used for minimization routines, it needs to calculate the negative of
    % the determinant.
end