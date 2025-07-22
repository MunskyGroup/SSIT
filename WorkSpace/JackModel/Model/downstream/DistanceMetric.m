function D = DistanceMetric(trueResults, testResults, TSBinIdx, distanceMetric)
    if nargin < 3 || isempty(TSBinIdx)
        TSBinIdx = [7];
    end
    
    if nargin < 4 || isempty(distanceMetric)
        distanceMetric = @JSDivergence;
    end

    allData = [trueResults.binCounts, testResults.binCounts];
    min_num_rna_in_any_bin = min(allData,[], 'all');
    max_num_rna_in_any_bin = max(allData,[], 'all');
    x = linspace(min_num_rna_in_any_bin, max_num_rna_in_any_bin);

    P = getPDF(trueResults);
    Q = getPDF(testResults);

    distances = zeros([length(TSBinIdx), size(P,2)]);
    for i = TSBinIdx % TS Bins
        for j = 1:size(P,2) % Spatial Bins
            p = P{i,j};
            q = Q{i,j};

            distances(i,j) = distanceMetric(p, q, x);
        end
    end
    
    D = sum(distances(:));
end