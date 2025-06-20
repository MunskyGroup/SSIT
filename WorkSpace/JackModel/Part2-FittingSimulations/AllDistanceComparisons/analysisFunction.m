function [distances, mean_distances, var_distances] = analysisFunction(results, name, parameter_range, row_index, titles_on, num_edges)
    num_eval_points = 100; % number of evaluation points for trapizoid method for calculating pdf overlap

    % Wasserstein
    % gets wasserstein distance
    distances = zeros(size(results,2), size(results,2), length(parameter_range), length(parameter_range));
    for j = 1:length(parameter_range)
        for k = 1:length(parameter_range)
            for n = 1:size(results,2)
                for n_1 = 1:size(results,2)
                    distances(n, n_1, j, k) = DistributionDistance_Wasserstein1D(results{j,n}, results{k,n_1}, num_edges);
                end
            end
        end
    end
    
    % gets mean_distances and vars from wasserstein
    mean_distances = zeros(length(parameter_range), length(parameter_range));
    var_distances = zeros(length(parameter_range), length(parameter_range));
    for j = 1:length(parameter_range)
        for k = 1:length(parameter_range)
            distances_ = distances(:,:,j,k);
            distances_ = distances_(~eye(size(distances_)));
            mean_distances(j,k) = mean(distances_);
            var_distances(j,k) = var(distances_);
        end
    end
    
    % get overlap between distance pdfs wasserstein
    overlap = zeros(length(parameter_range), length(parameter_range));
    for j = 1:length(parameter_range)
        self_distances = distances(:,:,j,j);
        self_distances = self_distances(~eye(size(self_distances)));
    
        for k = 1:length(parameter_range)
            deviant_distance = distances(:,:,j,k);
            deviant_distance = deviant_distance(~eye(size(deviant_distance)));
            overlap(j,k) = compute_overlap_area(deviant_distance, self_distances, num_eval_points);
        end
    end
    
    % Heatmap of mean wasterstein distance between parameters 
    figure;
    imagesc(log10(mean_distances));
    colorbar;
    xlabel(name);
    ylabel(name);
    if titles_on
        title(sprintf('Heatmap of log_{10} Distribution Distances for %s (Wasserstein1D)', name));
    end
    set(gca, 'XTick', 1:length(parameter_range), 'XTickLabel', parameter_range, ...
             'YTick', 1:length(parameter_range), 'YTickLabel', parameter_range);
    axis square;
    saveas(gcf, sprintf('%s_MeanDistancesHeatmap_Wasserstein1D.png', name));
    
    % Heatmap of Overlaps of PDF of distances between parameters
    figure;
    imagesc(overlap);
    colorbar;
    xlabel(name);
    ylabel(name);
    if titles_on
        title(sprintf('Heat map of overlap between Distance Distributions for Different %s (Wasserstein1D)', name));
    end
    set(gca, 'XTick', 1:length(parameter_range), 'XTickLabel', parameter_range, ...
             'YTick', 1:length(parameter_range), 'YTickLabel', parameter_range);
    axis square;
    saveas(gcf, sprintf('%s_Overlap_Wasserstein1D.png', name));
    
    std_across = sqrt(var_distances);
    upper = mean_distances + std_across;
    lower = mean_distances - std_across;
    
    % 1D Wasserstein Distance along row
    figure;
    hold on;
    xlabel(name)
    ylabel('Mean Distance ± STD');
    plot(parameter_range, mean_distances(row_index, :), 'LineWidth', 1.5, 'Color', [0 0 0.5])
    fill([parameter_range, fliplr(parameter_range)], [upper(row_index, :), fliplr(lower(row_index, :))], [0.9 0.9 1], 'EdgeColor', 'none');
    grid on;
    legend('±1 STD', 'Mean');
    saveas(gcf, sprintf('%s_MeanWithVariance_Wasserstein1D.png', name));
    

    % JSDivergence
    % get JSDivergence Distances 
    distances = zeros(size(results,2), size(results,2), length(parameter_range), length(parameter_range));
    for j = 1:length(parameter_range)
        for k = 1:length(parameter_range)
            for n = 1:size(results,2)
                for n_1 = 1:size(results,2)
                    distances(n, n_1, j, k) = DistributionDistance_JSDivergence(results{j,n}, results{k,n_1}, num_edges);
                end
            end
        end
    end
    
    % get mean and var of distances
    mean_distances = zeros(length(parameter_range), length(parameter_range));
    for j = 1:length(parameter_range)
        for k = 1:length(parameter_range)
            distances_ = distances(:,:,j,k);
            distances_ = distances_(~eye(size(distances_)));
            mean_distances(j,k) = mean(distances_);
            var_distances(j,k) = var(distances_);
        end
    end
    
    % get overlap between pdfs 
    overlap = zeros(length(parameter_range), length(parameter_range));
    for j = 1:length(parameter_range)
        self_distances = distances(:,:,j,j);
        self_distances = self_distances(~eye(size(self_distances)));
    
        for k = 1:length(parameter_range)
            deviant_distance = distances(:,:,j,k);
            deviant_distance = deviant_distance(~eye(size(deviant_distance)));
    
            overlap(j,k) = compute_overlap_area(deviant_distance, self_distances, num_eval_points);
        end
    end
    
    % Heatmap Mean Distance Heatmap JSDivergence
    figure;
    imagesc(log10(mean_distances));
    colorbar;
    xlabel(name);
    ylabel(name);
    set(gca, 'XTick', 1:length(parameter_range), 'XTickLabel', parameter_range, ...
             'YTick', 1:length(parameter_range), 'YTickLabel', parameter_range);
    axis square;
    saveas(gcf, sprintf('%s_MeanDistances_JSDivergence.png', name));
    
    % Heatmap Overlap JSDivergence
    figure;
    imagesc(overlap);
    colorbar;
    xlabel(name);
    ylabel(name);
    set(gca, 'XTick', 1:length(parameter_range), 'XTickLabel', parameter_range, ...
             'YTick', 1:length(parameter_range), 'YTickLabel', parameter_range);
    axis square;
    saveas(gcf, sprintf('%s_Overlap_JSDivergence.png', name));
    
    % 1D JSDivergence 
    std_across = sqrt(var_distances);
    upper = mean_distances + std_across;
    lower = mean_distances - std_across;
    
    figure;
    hold on;
    xlabel(name)
    ylabel('Mean Distance ± STD');
    plot(parameter_range, mean_distances(row_index, :), 'LineWidth', 1.5, 'Color', [0 0 0.5])
    fill([parameter_range, fliplr(parameter_range)], [upper(row_index, :), fliplr(lower(row_index, :))], [0.9 0.9 1], 'EdgeColor', 'none');
    grid on;
    legend('±1 STD', 'Mean');
    saveas(gcf, sprintf('%s_MeanWithVariance_JSDivergence.png', name));

end


function P = Prob(results, edges)
    % Calculates the probability density function as matrix with edges
    % as bins
    % the rows of the matrix corrisopond to 
    % 1-5 - cells that the TS was identifiable 
    % 6 - cells with uknown TS
    % 7 - all cells

    % collumns of along x
    % bins from left to right
    if nargin < 2 || isempty(edges)
        edges = [0,10,20,30,40,50,60,70,80,90,100];
    end

    RNACountBins = results.binCounts; 
    TSbins = results.binTS;

    P = cell([7, size(RNACountBins, 1)]);

    % get probabilites 
    for ibinloc = 1:5
        iCellsInTSBin = (TSbins == ibinloc);
        for iBin = 1:size(RNACountBins, 1)
            data = RNACountBins(iBin, iCellsInTSBin);
            if isempty(data)
                pdf = zeros(1, length(edges)-1);
            else
                pdf = histcounts(data, 'BinEdges', edges, 'Normalization', 'probability');
            end
            P{ibinloc, iBin} = pdf; % pdf of bined TS 
        end
    end

    iCellsInTSBin = (TSbins == 0);
    for iBin = 1:size(RNACountBins, 1)
        data = RNACountBins(iBin, iCellsInTSBin);
        if isempty(data)
            pdf = zeros(1, length(edges)-1);
        else
            pdf = histcounts(data, 'BinEdges', edges, 'Normalization', 'probability');
        end
        P{6, iBin} = pdf;  % pdf of unknown ts sites
    end

    for iBin = 1:size(RNACountBins, 1)
        data = RNACountBins(iBin, :);
        if isempty(data)
            pdf = zeros(1, length(edges)-1);
        else
            pdf = histcounts(data, 'BinEdges', edges, 'Normalization', 'probability');
        end
        P{7, iBin} = pdf; % pdf of all cells regardless of TS
    end
end


function d = Wasserstein1D(P, Q, binCenters)
    % 1D Wasserstein (Earth Mover's) distance
    cdfP = cumsum(P);
    cdfQ = cumsum(Q);
    d = sum(abs(cdfP - cdfQ) .* diff([0 binCenters]));
end


function d = JSDivergence(P, Q)
    epsilon = 1e-12;
    P = P + epsilon;
    Q = Q + epsilon;

    P = P / sum(P);
    Q = Q / sum(Q);

    M = 0.5 * (P + Q);
    d = 0.5 * sum(P .* log2(P ./ M)) + 0.5 * sum(Q .* log2(Q ./ M));
end


function metric = DistributionDistance_Wasserstein1D(results1, results2, num_bins)
    % auto calculates edges 
    allData = [results1.binCounts, results2.binCounts];
    min_num_rna_in_any_bin = min(allData,[], 'all');
    max_num_rna_in_any_bin = max(allData,[], 'all');
    edges = linspace(min_num_rna_in_any_bin, max_num_rna_in_any_bin, num_bins);
    binCenters = edges(1:end-1) + diff(edges)/2;

    if max_num_rna_in_any_bin ~= 0
        % calculate pdfs vector with same edges 
        P = Prob(results1, edges);
        Q = Prob(results2, edges);
    
        % calculate metric
        distances = zeros(size(P));
        for i = 1:size(P,1)
            for j = 1:size(P,2)
                p = P{i,j};
                q = Q{i,j};
    
                distances(i,j) = Wasserstein1D(p, q, binCenters);
            end
        end
        
        metric =  sum(distances(7,:));
    else
        metric = 0;
    end
end

function metric = DistributionDistance_JSDivergence(results1, results2, num_bins)
    % auto calculates edges 
    allData = [results1.binCounts, results2.binCounts];
    min_num_rna_in_any_bin = min(allData,[], 'all');
    max_num_rna_in_any_bin = max(allData,[], 'all');
    edges = linspace(min_num_rna_in_any_bin, max_num_rna_in_any_bin, num_bins);

    if max_num_rna_in_any_bin ~= 0
        % calculate pdfs vector with same edges 
        P = Prob(results1, edges);
        Q = Prob(results2, edges);
    
        % calculate metric
        distances = zeros(size(P));
        for i = 1:size(P,1)
            for j = 1:size(P,2)
                p = P{i,j};
                q = Q{i,j};
    
                % p = p / sum(p + eps);
                % q = q / sum(q + eps);
    
                distances(i,j) = JSDivergence(p, q);
            end
        end
    
        metric = sum(distances(7,:));
    else
        metric = 0;
    end
end

function [overlap_area] = compute_overlap_area(samples1, samples2, num_eval_points)
    % Step 1: Define common x range
    x_min = min([samples1; samples2]);
    x_max = max([samples1; samples2]);
    x = linspace(x_min, x_max, num_eval_points);  % Common evaluation points

    % Step 2: Estimate PDFs using KDE
    [pdf1, x1] = ksdensity(samples1, x);
    [pdf2, x2] = ksdensity(samples2, x);

    % Step 3: Compute pointwise minimum
    min_pdf = min(pdf1, pdf2);

    % Step 4: Integrate overlap area
    overlap_area = trapz(x, min_pdf);
end