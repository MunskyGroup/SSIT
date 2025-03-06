classdef FIMMetricKind
    enumeration
        Determinant % Maximize the expected determinant of the FIM
        DeterminantCoveriance % Minimize the expected determinant of MLE covariance 
        SmallestEigenvalue % Maximize the smallest eigenvalue of the FIM
        Trace % Maximize the trace of the FIM
    end
end