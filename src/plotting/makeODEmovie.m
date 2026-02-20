function makeODEmovie(ODE_soln, speciesNames, timeVec, videoFileName)
    % plotODE_video - Animates and saves ODE trajectories as a video.
    %
    % Inputs:
    %   * ODE_soln - struct with field 'ode' (nTime × nSpecies)
    %   * speciesNames (optional) - cell array of species names
    %   * timeVec (optional) - time vector (nTime × 1)
    %   * videoFileName (optional) - output video name (default: 'ode_trajectories.mp4')

    X = ODE_soln.ode; % [nTime × nSpecies]
    [nTime, numSpecies] = size(X);

    if nargin < 4 || isempty(videoFileName)
        videoFileName = 'ode_trajectories.mp4';
    end
    if nargin < 3 || isempty(timeVec)
        timeVec = 1:nTime;
    end
    if nargin < 2 || isempty(speciesNames)
        speciesNames = arrayfun(@(s) sprintf('Species %d', s), 1:numSpecies, 'UniformOutput', false);
    elseif length(speciesNames) ~= numSpecies
        error('The number of species names must match the number of species (%d).', numSpecies);
    end

    % Prepare video writer
    v = VideoWriter(videoFileName, 'MPEG-4');
    v.FrameRate = 10;
    open(v);

    % Setup figure
    figure;
    hold on;
    colors = lines(numSpecies);

    h = gobjects(1, numSpecies);
    for s = 1:numSpecies
        h(s) = plot(NaN, NaN, '-', 'LineWidth', 2, 'Color', colors(s, :));
    end

    xlabel('Time');
    ylabel('Molecule Count / Concentration');
    title('ODE Trajectories Animation');
    legend(speciesNames, 'Location', 'southeast');
    grid on;

    % Animate each time point
    for tIdx = 2:nTime
        tNow = timeVec(1:tIdx);
        for s = 1:numSpecies
            set(h(s), 'XData', tNow, 'YData', X(1:tIdx, s));
        end
        drawnow;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    close(v);
    disp(['Video saved to ', videoFileName]);
end
