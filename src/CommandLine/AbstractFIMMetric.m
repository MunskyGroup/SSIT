classdef AbstractFIMMetric
    %AbstractFIMMetric calculates an individual function (e.g., trace or
    %determinant) of an FIM upon demand.

    methods (Access = protected)
        function submatrix = getSubmatrixOfFIM(obj, FIM)
            % By default, all entries within the FIM will be considered, so
            % the submatrix equals the full matrix. Descendants can
            % override.

            submatrix = FIM;
        end
    end % Protected methods

    methods (Abstract, Access = protected)
        metricResult = calculateForFIMInternal(obj, FIM)
    end % Abstract protected methods

    methods
        function metricResult = calculateForFIM(obj, FIM)
            submatrix = getSubmatrixOfFIM(obj, FIM);
            metricResult = obj.calculateForFIMInternal(submatrix);
        end
    end % Public methods
end