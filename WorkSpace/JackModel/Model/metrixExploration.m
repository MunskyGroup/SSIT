clear 
close all

%% Data to be analysised
file_location = ["C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00325_kr-86250_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00325_kr-75000_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00325_kr-63750_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00325_kr-52500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00325_kr-97500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002875_kr-97500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002875_kr-86250_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002875_kr-75000_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002875_kr-63750_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002875_kr-52500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.0025_kr-97500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.0025_kr-86250_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.0025_kr-75000_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.0025_kr-63750_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.0025_kr-52500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002125_kr-97500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002125_kr-86250_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002125_kr-75000_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002125_kr-63750_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.002125_kr-52500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00175_kr-97500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00175_kr-86250_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00175_kr-75000_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00175_kr-63750_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"
"C:\Users\Jack\Documents\GitHub\HeLa_Nuclear_RNA_Modeling\codes\matlab\RNA_Modeling\simulation_kon-0.0007_koff-5.25e-05_w-0.00175_kr-52500_D1-0.01_D2-5_D3-4_gam1-0.035_gam2-0.0025_gam3-0.001.csv"];


%% Ripley and pair correlation
for i = 1:length(file_location)
    T = readtable(file_location(i));
    for c = 1:length(unique(T.cell_id))
    cell = T(T.cell_id == c, :);
    x = table2array(cell(:, 'x'));
    y = table2array(cell(:, 'y'));

    % Parameters
    r = 1:1:20; % Distance range for K-function
    binEdges = r; % Bin edges for histogram
    
    % Compute Ripley's K-function
    K = ripleyK(x, y, r, binEdges);
    
    % Plot K-function
    figure;
    plot(r, K, '-o');
    xlabel('Distance (r)');
    ylabel('K(r) - r');
    title(sprintf('Ripleys K-function %s', file_location(i)));
    
    % Compute pair correlation function
    g = pairCorrelation(x, y, r, binEdges);
    
    % Plot pair correlation function
    figure;
    plot(r, g, '-o');
    xlabel('Distance (r)');
    ylabel('g(r)');
    title(sprintf('Pair Correlation Function %s', file_location(i)));

    end

end





%% Metrics 
function K = ripleyK(x, y, r, binEdges)
    % x, y: Coordinates of points
    % r: Maximum distance to consider
    % binEdges: Bin edges for distance intervals

    numPoints = length(x);
    distances = pdist2([x, y], [x, y]);
    
    % Remove self-distances (diagonal elements)
    distances = distances + eye(numPoints) * Inf;
    
    % Compute the count of distances within each bin
    counts = histcounts(distances(:), binEdges);
    
    % Compute the area of annular rings for each bin
    area = pi * (binEdges(2:end).^2 - binEdges(1:end-1).^2); % Area difference for each bin
    density = numPoints / (pi * r(end)^2); % Density of points
    expected = density * area * numPoints;
    
    % Compute Ripley's K-function
    K = (counts ./ expected) - 1;
end


function g = pairCorrelation(x, y, r, binEdges)
    % x, y: Coordinates of points
    % r: Maximum distance to consider
    % binEdges: Bin edges for distance intervals

    numPoints = length(x);
    distances = pdist2([x, y], [x, y]);
    
    % Remove self-distances (diagonal elements)
    distances = distances + eye(numPoints) * Inf;
    
    % Compute the count of distances within each bin
    counts = histcounts(distances(:), binEdges);
    
    % Compute the area of annular rings for each bin
    area = pi * (binEdges(2:end).^2 - binEdges(1:end-1).^2); % Area difference for each bin
    density = numPoints / (pi * r(end)^2); % Density of points
    expected = density * area * numPoints;
    
    % Compute pair correlation function
    g = counts ./ expected;
end
