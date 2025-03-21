function plotODE(ODE_GR_soln, speciesNames, timeVec)
    % plotODE - Plots ODE solution for all species over time.
    %
    % Inputs:
    %   ODE_GR_soln: Struct with field 'ode' (nTime × nSpecies)
    %   speciesNames (optional): Cell array of species names
    %   timeVec (optional): Time vector [nTime × 1]

    X = ODE_GR_soln.ode;  % size: [nTime × nSpecies]
    [nTime, numSpecies] = size(X);

    % Default time vector if not provided
    if nargin < 3 || isempty(timeVec)
        timeVec = 1:nTime;  % Use indices if no time vector
    end

    % Generate default species names if not provided
    if nargin < 2 || isempty(speciesNames)
        speciesNames = arrayfun(@(s) sprintf('Species %d', s), 1:numSpecies, 'UniformOutput', false);
    elseif length(speciesNames) ~= numSpecies
        error('The number of species names must match the number of species (%d).', numSpecies);
    end

    % Plot
    figure; hold on;
    colors = lines(numSpecies);
    for s = 1:numSpecies
        plot(timeVec, X(:, s), 'LineWidth', 2, 'Color', colors(s, :));
    end

    xlabel('Time');
    ylabel('Molecule Count / Concentration');
    title('ODE Solution Trajectories');
    legend(speciesNames, 'Location', 'Best');
    grid on;
    hold off;
end