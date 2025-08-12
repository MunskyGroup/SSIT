function plotODE(ODE_soln, speciesNames, timeVec)
    % plotODE - Plots ODE solution for all model species over time.
    %
    % Inputs:
    %   * ODE_soln - struct with field 'ode' (nTime × nSpecies)
    %   * speciesNames (optional) - cell array of species names for plot
    %                               legend
    %   * timeVec (optional) - time vector [nTime × 1]
    %
    % Example: plotODE(Model_ODEsoln,Model_ODE.species,Model_ODE.tSpan)

    X = ODE_soln.ode;  % size: [nTime × nSpecies]
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
        plot(timeVec, X(:, s), 'LineWidth', 4, 'Color', colors(s, :));
    end

    xlabel('Time');
    ylabel('Molecule Count / Concentration');
    title('ODE Solution Trajectories');
    legend(speciesNames, 'Location', 'Best');
    grid on;
    hold off;
end