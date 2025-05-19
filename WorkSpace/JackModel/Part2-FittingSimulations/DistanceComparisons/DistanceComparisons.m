%% Test distance meterics in several directions
clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'))

s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

% Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot)


%% Parameter Set 1
makePlot = false;

N=10;

NCells = 2000;
kon = 1e-3;
koff = 7.5e-5;
w = 0.0025;
kex = 750;
kr = 7.5e4;
D = [0.01,5,4];
gam =[0.035;0.0025;0.001];

if false
    results = cell(N);
    for n = 1:N
        results{n} = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot);
    end
    save("PartOne", "results")
end

results = load("PartOne.mat", "results").results;

n = length(results);
distances = zeros(n);

for i = 1:n
    for j = 1:n
        p = results{i};
        q = results{j};
        distances(i, j) = DistributionDistance(p, q);
    end
end

self_distances = distances(~eye(size(distances)));

%% Parameter Set 2
kon = 1e-4; %% 
koff = 7.5e-5;
w = 0.0025;
kex = 750;
kr = 7.5e4;
D = [0.01,5,4];
gam =[0.035;0.0025;0.001];

if false
    results1 = cell(N);
    for n = 1:N
        results1{n} = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot);
    end
    save("PartTwo", "results1")
end

results1 = load("PartTwo.mat", "results1").results1;

n = length(results1);
distances = zeros(n);
distances1 = zeros(n);

for i = 1:n
    for j = 1:n
        p = results1{i};
        q = results1{j};
        distances(i, j) = DistributionDistance(p, q);

        p = results{i};
        q = results1{j};
        distances1(i, j) = DistributionDistance(p, q);
    end
end

self_distances1 = distances(~eye(size(distances)));
slightly_different_distances1 = distances1;

%% Parameter Set 3
kon = 1; %% 
koff = 7.5e-5;
w = 0.0025;
kex = 750;
kr = 7.5e4;
D = [0.01,5,4];
gam =[0.035;0.0025;0.001];

if false
    results2 = cell(N);
    for n = 1:N
        results2{n} = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot);
    end
    save("PartThree", "results2")
end

results2 = load("PartThree.mat", "results2").results2;

n = length(results2);
distances = zeros(n);
distances1 = zeros(n);
distances2 = zeros(n);

for i = 1:n
    for j = 1:n
        p = results2{i};
        q = results2{j};
        distances(i, j) = DistributionDistance(p, q);

        p = results{i};
        q = results2{j};
        distances1(i, j) = DistributionDistance(p, q);

        p = results1{i};
        q = results2{j};
        distances2(i, j) = DistributionDistance(p, q);
    end
end

self_distances2 = distances(~eye(size(distances)));
very_different_distances2 = distances1;
very_different_distances3 = distances2;




%% Display First Parts 1-3
figure;
hold on;

histogram(self_distances, 'BinMethod', 'auto', ...
    'Normalization', 'probability', ...
    'FaceColor', [0.60 0.80 0.35], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % light green

histogram(self_distances1, 'BinMethod', 'auto', ...
    'Normalization', 'probability', ...
    'FaceColor', [0.40 0.70 0.20], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % medium green

histogram(self_distances2, 'BinMethod', 'auto', ...
    'Normalization', 'probability', ...
    'FaceColor', [0.20 0.50 0.10], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % dark green

% Between-group comparisons (distinct colors)
histogram(slightly_different_distances1, 'BinMethod', 'auto', ...
    'Normalization', 'probability', ...
    'FaceColor', [0.49 0.18 0.56], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % purple

histogram(very_different_distances2, 'BinMethod', 'auto', ...
    'Normalization', 'probability', ...
    'FaceColor', [0.30 0.75 0.93], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % light blue

histogram(very_different_distances3, 'BinMethod', 'auto', ...
    'Normalization', 'probability', ...
    'FaceColor', [0.85 0.33 0.10], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % reddish-orange

% Set log scale for x-axis
set(gca, 'XScale', 'log');

% Title and axis labels
title('Wasserstein Distance Comparisons');
xlabel('Wasserstein (Log Scale)');
ylabel('Probability');

% Updated legend
legend('Distance(ParamSet1, ParamSet1)', ...
       'Distance(ParamSet2, ParamSet2)', ...
       'Distance(ParamSet3, ParamSet3)', ...
       'Distance(ParamSet1, ParamSet2)', ...
       'Distance(ParamSet2, ParamSet3)', ...
       'Distance(ParamSet1, ParamSet3)', ...
       'Location', 'best');

hold off;


%% Distance Functions
function P = Prob(results, edges)
    if nargin < 2 || isempty(edges)
        edges = [0,10,20,30,40,50,60,70,80,90,100];
    end

    RNACountBins = results.binCounts; 
    TSbins = results.binTS;

    P = cell([7, size(RNACountBins, 1)]);

    for ibinloc = 1:5
        iCellsInTSBin = (TSbins == ibinloc);
        for iBin = 1:size(RNACountBins, 1)
            data = RNACountBins(iBin, iCellsInTSBin);
            if isempty(data)
                pdf = zeros(1, length(edges)-1);
            else
                pdf = histcounts(data, 'BinEdges', edges, 'Normalization', 'probability');
            end
            P{ibinloc, iBin} = pdf;
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
        P{6, iBin} = pdf;
    end

    for iBin = 1:size(RNACountBins, 1)
        data = RNACountBins(iBin, :);
        if isempty(data)
            pdf = zeros(1, length(edges)-1);
        else
            pdf = histcounts(data, 'BinEdges', edges, 'Normalization', 'probability');
        end
        P{7, iBin} = pdf;
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

function d = Wasserstein1D(P, Q, binCenters)
    % 1D Wasserstein (Earth Mover's) distance
    cdfP = cumsum(P);
    cdfQ = cumsum(Q);
    d = sum(abs(cdfP - cdfQ) .* diff([0 binCenters]));
end


function metric = DistributionDistance(results1, results2, edges)

    if nargin < 3 || isempty(edges)
        allData = [results1.binCounts(:); results2.binCounts(:)];
        edges = linspace(min(allData), max(allData), 21);
    end

    binCenters = edges(1:end-1) + diff(edges)/2;

    P = Prob(results1, edges);
    Q = Prob(results2, edges);

    % Compare overall RNA count distribution
    data1 = sum(results1.binCounts, 1);
    data2 = sum(results2.binCounts, 1);

    pdf1 = histcounts(data1, edges, 'Normalization', 'probability');
    pdf2 = histcounts(data2, edges, 'Normalization', 'probability');

    countDistance = Wasserstein1D(pdf1, pdf2, binCenters);

    distances = zeros(size(P));
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            p = P{i,j};
            q = Q{i,j};

            % Pad if necessary
            if isempty(p); p = zeros(1, length(binCenters)); end
            if isempty(q); q = zeros(1, length(binCenters)); end

            p = p / sum(p + eps);
            q = q / sum(q + eps);

            distances(i,j) = Wasserstein1D(p, q, binCenters);
        end
    end

    metric = countDistance + sum(distances(7,:));
end



