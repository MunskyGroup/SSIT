function [] = plotMarginalInSensTab(app)
% this function plots the marginals generated from the Sensitivity Tab in a
% new popout window

% Find the time index to plot
T_array = eval(app.SensPrintTimesEditField.Value);
[~,j] = min(abs(T_array-app.SensPlotTimeSlider.Value));

%% Compute the marginal distributions
solutionFormat = app.SensFspTabOutputs.solutions.format;

if (strcmp(solutionFormat, 'forward'))    
    mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states, ...
                              app.SensFspTabOutputs.solutions.data{j}.p,false);
%     sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states,...
%         app.SensFspTabOutputs.solutions.data{j}.dp(:, ipar),false);
else    
    mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.ps{j}.states, ...
                              app.SensFspTabOutputs.solutions.data.ps{j}.p,false);
%     sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.dps{j, ipar}.states,...
%         app.SensFspTabOutputs.solutions.data.dps{j, ipar}.p,false);
end

%% Plot the marginal distribution
speciesToPlot = [];
legends = {};
for iSpecies = 1:3
    plotThisSpecies = eval(['app.SensMarginalX', num2str(iSpecies), 'CheckBox.Value']);
    if plotThisSpecies
        speciesToPlot = [speciesToPlot iSpecies];
        legends = [legends char(app.NameTable.Data(iSpecies,2))];
    end
    
end

%% Plot trajectories, each species have a color
speciesColors = {'b', 'r', 'g'};
for iSim = 1:3
    % Plot the trajectories
    for iSim = speciesToPlot
        stairs([0:length(mdist{iSim})], [mdist{iSim};0], speciesColors{iSim},'Linewidth',2);
        hold on
    end
end
%% Add Legend Title and Axes
legend(legends,'Location','best');
title(sprintf('Marginals at time: t = %1.2f',T_array(j)));
xlabel('Species Count');
ylabel('Probability');
set(gca,'Fontsize',20);
end


