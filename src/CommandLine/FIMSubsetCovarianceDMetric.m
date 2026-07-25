classdef FIMSubsetCovarianceDMetric < FIMCovarianceDMetric & ...
        AbstractFIMSubsetMetric
    
    properties
        NumberOfTimes (1, 1) double {mustBeInteger, mustBePositive} = 1
    end

    methods (Access = protected)
        function submatrix = getSubmatrixOfFIM(obj, FIM)
            % This override returns the subset of the FIM corresponding to
            % the specified parameter indices. Descendants should therefore
            % expect that this is the full FIM for purposes of metric
            % calculation, unless further or alternative overrides are
            % implemented.

            k = obj.ParameterIndices;
            ek = zeros(length(k), obj.NumberOfTimes);
            ek(1:length(k), k) = eye(length(k));

            submatrix = ek * inv(FIM) * ek';
        end
    end % Protected methods
end