classdef AbstractFIMSubsetMetric < AbstractFIMMetric
    properties
        ParameterIndices (1, :) {mustBeInteger, mustBePositive} = []
    end

    methods (Access = protected)
        function submatrix = getSubmatrixOfFIM(obj, FIM)
            % This override returns the subset of the FIM corresponding to
            % the specified parameter indices. Descendants should therefore
            % expect that this is the full FIM for purposes of metric
            % calculation, unless further or alternative overrides are
            % implemented.

            submatrix = FIM(obj.ParameterIndices, obj.ParameterIndices);
        end
    end % Protected methods

    methods
        function obj = set.ParameterIndices(obj, indices)
            obj.ParameterIndices = unique(sort(indices));
        end
    end % Public methods
end