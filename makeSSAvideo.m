function makeSSAvideo(ssaSoln, speciesIdx, numTraj, speciesNames, videoFileName)
    % plotSSA_video - Creates a video of SSA trajectories over time.
    %
    % Inputs:
    %   * ssaSoln - struct containing SSA simulation results
    %   * speciesIdx - index of the species to plot (1 to N) or 'all'
    %   * numTraj - number of trajectories to display
    %   * speciesNames (optional) - cell array of species names
    %   * videoFileName (optional) - name of the .mp4 file to save (default: 'ssa_trajectories.mp4')

    if nargin < 5 || isempty(videoFileName)
        videoFileName = 'ssa_trajectories.mp4';
    end

    numSpecies = size(ssaSoln.trajs, 1);
    numTotalTraj = size(ssaSoln.trajs, 3);
    numTraj = min(numTraj, numTotalTraj);

    if nargin < 4 || isempty(speciesNames)
        speciesNames = arrayfun(@(s) sprintf('Species %d', s), 1:numSpecies, 'UniformOutput', false);
    elseif length(speciesNames) ~= numSpecies
        error('The number of species names must match the number of species (%d).', numSpecies);
    end

    T = ssaSoln.T_array;
    validIdx = T >= 0;
    T = T(validIdx);
    trajs = ssaSoln.trajs(:, validIdx, :);

    speciesColors = lines(numSpecies);
    randIdx = randperm(numTotalTraj, numTraj); % pick trajectories to show

    % Setup video writer
    v = VideoWriter(videoFileName, 'MPEG-4');
    v.FrameRate = 10;
    open(v);

    % Create figure for plotting
    figure;
    hold on;

    % Initialize plot handles
    if strcmp(speciesIdx, 'all')
        h = gobjects(numSpecies, numTraj);
        for s = 1:numSpecies
            for i = 1:numTraj
                h(s, i) = plot(NaN, NaN, '-', 'Color', [speciesColors(s, :) 0.3]);
            end
        end
        meanLines = gobjects(1, numSpecies);
        for s = 1:numSpecies
            meanLines(s) = plot(NaN, NaN, 'Color', speciesColors(s, :), 'LineWidth', 2);
        end
    else
        s = speciesIdx;
        h = gobjects(1, numTraj);
        for i = 1:numTraj
            h(i) = plot(NaN, NaN, '-', 'Color', [speciesColors(s, :) 0.3]);
        end
        meanLine = plot(NaN, NaN, 'Color', speciesColors(s, :), 'LineWidth', 2);
    end

    xlabel('Time');
    ylabel('Molecule Count');
    if strcmp(speciesIdx, 'all')
        legend(meanLines, speciesNames, 'Location', 'Best');
    else
        legend(meanLine, speciesNames{speciesIdx}, 'Location', 'Best');
    end
    grid on;

    % Animate over time
    for tIdx = 2:length(T)
        tNow = T(1:tIdx);

        if strcmp(speciesIdx, 'all')
            for s = 1:numSpecies
                Xs = squeeze(trajs(s, 1:tIdx, randIdx));
                for i = 1:numTraj
                    set(h(s, i), 'XData', tNow, 'YData', Xs(:, i));
                end
                set(meanLines(s), 'XData', tNow, 'YData', mean(Xs, 2));
            end
        else
            Xs = squeeze(trajs(speciesIdx, 1:tIdx, randIdx));
            for i = 1:numTraj
                set(h(i), 'XData', tNow, 'YData', Xs(:, i));
            end
            set(meanLine, 'XData', tNow, 'YData', mean(Xs, 2));
        end

        drawnow;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    close(v);
    disp(['Video saved to ', videoFileName]);
end