classdef MBExperimentFIMOptimizedDesignStrategyInfo < ...
        MBExperimentDesignStrategyInfo
    %MBExperimentFIMOptimizedDesignStrategyInfo contains
    %information needed by all FIM-optimized design strategies to design an
    %experiment.

    properties
        CovariancePrior (1, :) double = []
        FIMs (:, :) cell
        %   * 'Nc' - an optimal guess for the optimal experiment
        %            design
        Metric (1, 1) AbstractFIMMetric = FIMDMetric();
        NcGuess = []
        %   * 'NcFixed' - a minimal number of cells to measure at each
        %      time point; this is useful for subsequent experiment
        %      design, having already obtained measured cells from a
        %      previous experiment
        NcFixed = []
        %   * 'NcMax' - maximum total number of cells allowed for each
        %      time point; this is useful in simulated experiment design,
        %      where there are only so many cells available in the real
        %      data
        NcMax = []
        Statistic (1, 1) FIMMetricStatistic = ...
            FIMMetricStatistic.Mean;
    end
end