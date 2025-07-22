function P = getPDF(results)
% Calculates the probability density function as matrix with edges
% as bins
% the rows of the matrix corrisopond to
% 1-5 - cells that the TS was identifiable
% 6 - cells with uknown TS
% 7 - all cells

RNACountBins = results.binCounts;
TSbins = results.binTS;

P = cell([7, size(RNACountBins, 1)]);

% Main loop for TS bins 1-5
for ibinloc = 1:5
    iCellsInTSBin = (TSbins == ibinloc);
    for iBin = 1:size(RNACountBins, 1)
        data = RNACountBins(iBin, iCellsInTSBin);

        [pdf,edges] = histcounts(data, 'Normalization', 'probability');
        binCenters = (edges(1:end-1) + edges(2:end)) / 2;

        if numel(binCenters) > 1
            interpFcn = griddedInterpolant(binCenters, pdf);
        else
            interpFcn = @(x) zeros(size(x));  % Handle case with insufficient data
        end

        P{ibinloc, iBin} = @(x) max(interpFcn(x), 0);  % ensure non-negative output
    end
end

% Unknown TS bin (TSbin == 0)
iCellsInTSBin = (TSbins == 0);
for iBin = 1:size(RNACountBins, 1)
    data = RNACountBins(iBin, iCellsInTSBin);

    [pdf,edges] = histcounts(data, 'Normalization', 'probability');
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;

    if numel(binCenters) > 1
        interpFcn = griddedInterpolant(binCenters, pdf);
    else
        interpFcn = @(x) zeros(size(x));  % Handle case with insufficient data
    end

    P{6, iBin} = @(x) max(interpFcn(x), 0);
end

% All cells
for iBin = 1:size(RNACountBins, 1)
    data = RNACountBins(iBin, :);

    [pdf,edges] = histcounts(data, 'Normalization', 'probability');
    binCenters = (edges(1:end-1) + edges(2:end)) / 2;

    if numel(binCenters) > 1
        interpFcn = griddedInterpolant(binCenters, pdf);
    else
        interpFcn = @(x) zeros(size(x));  % Handle case with insufficient data
    end

    P{7, iBin} = @(x) max(interpFcn(x), eps);
end
end