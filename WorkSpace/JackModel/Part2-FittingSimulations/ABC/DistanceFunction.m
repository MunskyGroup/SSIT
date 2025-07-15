function metric = DistanceFunction(results1, results2, num_bins)
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


function d = JSDivergence(P, Q)
    epsilon = 1e-12;
    P = P + epsilon;
    Q = Q + epsilon;

    P = P / sum(P);
    Q = Q / sum(Q);

    M = 0.5 * (P + Q);
    d = 0.5 * sum(P .* log2(P ./ M)) + 0.5 * sum(Q .* log2(Q ./ M));
end




