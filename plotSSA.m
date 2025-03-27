function plotSSA(ssaSoln, speciesIdx, numTraj, speciesNames)
    % plotSSA - Plots SSA trajectories and histograms from ssaSoln struct.
    %
    % Inputs:
    %   ssaSoln: Struct containing SSA simulation results
    %   speciesIdx: Index of the species to plot (1 to N) or 'all' to plot all species
    %   numTraj: Number of trajectories to display (max available)
    %   speciesNames (optional): Cell array of species names (must match numSpecies)

    numSpecies = size(ssaSoln.trajs, 1); % Automatically detect number of species
    numTotalTraj = size(ssaSoln.trajs, 3);
    numTraj = min(numTraj, numTotalTraj); % Ensure we don't exceed available trajectories

    % If species names are not provided, generate default names
    if nargin < 4 || isempty(speciesNames)
        speciesNames = arrayfun(@(s) sprintf('Species %d', s), 1:numSpecies, 'UniformOutput', false);
    elseif length(speciesNames) ~= numSpecies
        error('The number of species names must match the number of species (%d).', numSpecies);
    end

    % Extract time points and filter valid ones (t >= 0)
    T = ssaSoln.T_array;
    validIdx = T >= 0;
    T = T(validIdx);
    
    % Locate the index closest to time = 100
    [~, t100_idx] = min(abs(T - 100));

    % Define colors dynamically based on the number of species
    speciesColors = lines(numSpecies); % Generate distinct colors for all species
    
    figure; hold on;
    legendEntries = {}; % Store legend labels
    legendHandles = []; % Store handles for legend colors

    if strcmp(speciesIdx, 'all')
        % Plot all species in different colors
        for s = 1:numSpecies
            X = squeeze(ssaSoln.trajs(s, validIdx, :)); % Extract valid species trajectories
            randIdx = randperm(numTotalTraj, numTraj); % Select random trajectories
            
            for i = 1:numTraj
                plot(T, X(:, randIdx(i)), 'Color', [speciesColors(s, :), 0.2]); % Transparent individual trajectories
            end
            
            % Plot mean trajectory in the correct color
            h = plot(T, mean(X, 2), 'Color', speciesColors(s, :), 'LineWidth', 2);
            
            % Store handle and label for legend
            legendHandles(end+1) = h; %#ok<AGROW>
            legendEntries{end+1} = speciesNames{s}; % Use provided species name
        end
    else
        % Single species case
        if speciesIdx < 1 || speciesIdx > numSpecies
            error('speciesIdx must be between 1 and %d, or ''all''.', numSpecies);
        end
        
        X = squeeze(ssaSoln.trajs(speciesIdx, validIdx, :)); % Extract valid species trajectories
        randIdx = randperm(numTotalTraj, numTraj); % Select random trajectories
        
        for i = 1:numTraj
            plot(T, X(:, randIdx(i)), 'Color', [speciesColors(speciesIdx, :), 0.2]); % Transparent individual trajectories
        end
        
        % Plot mean trajectory in the correct color
        h = plot(T, mean(X, 2), 'Color', speciesColors(speciesIdx, :), 'LineWidth', 2);
        
        % Store handle and label for legend
        legendHandles = h;
        legendEntries = {speciesNames{speciesIdx}};
    end

    % Labels
    xlabel('Time');
    ylabel('Molecule Count');
    if strcmp(speciesIdx, 'all')
        title('SSA Trajectories for All Species (Starting at t=0)');
    else
        title(sprintf('SSA Trajectories for %s (Starting at t=0)', speciesNames{speciesIdx}));
    end
    legend(legendHandles, legendEntries, 'Location', 'Best'); % Ensure correct species names in legend
    grid on;
    hold off;

    % -------- HISTOGRAM AT TIME = 100 --------
    figure;
    numRows = ceil(sqrt(numSpecies)); % Adjust subplot grid dynamically
    numCols = ceil(numSpecies / numRows);

    if strcmp(speciesIdx, 'all')
        for s = 1:numSpecies
            subplot(numRows, numCols, s);
            X_t100 = squeeze(ssaSoln.trajs(s, t100_idx, :)); % Extract values at t=100
            histogram(X_t100, 'FaceColor', speciesColors(s, :), 'EdgeColor', 'k');
            xlabel(sprintf('Molecule Count (%s)', speciesNames{s}));
            ylabel('Frequency');
            title(sprintf('Distribution at Time = 100 (%s)', speciesNames{s}));
            grid on;
        end
    else
        X_t100 = squeeze(ssaSoln.trajs(speciesIdx, t100_idx, :)); % Extract values at t=100
        histogram(X_t100, 'FaceColor', speciesColors(speciesIdx, :), 'EdgeColor', 'k');
        xlabel(sprintf('Molecule Count (%s)', speciesNames{speciesIdx}));
        ylabel('Frequency');
        title(sprintf('Distribution at Time = 100 (%s)', speciesNames{speciesIdx}));
        grid on;
    end
end