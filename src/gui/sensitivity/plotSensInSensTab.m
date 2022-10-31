function [] = plotSensInSensTab(app)
% this function plots the Sensitivity generated from the Sensitivity Tab in a
% new popout window

% Find the time index to plot
T_array = eval(app.SensPrintTimesEditField.Value);
[~,j] = min(abs(T_array-app.SensPlotTimeSlider.Value));

% Find the parameter index to plot
ipar = 1;
num_pars = length(app.ReactionsTabOutputs.parameters);
while (~strcmp(app.ReactionsTabOutputs.parameters{ipar}, app.SensParDropDown.Value) && ipar < num_pars)
    ipar = ipar + 1;
end

%% Compute the marginal distributions
solutionFormat = app.SensFspTabOutputs.solutions.format;

if (strcmp(solutionFormat, 'forward'))    
    sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states,...
        app.SensFspTabOutputs.solutions.data{j}.dp(:, ipar),false);
else    
    sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.dps{j, ipar}.states,...
        app.SensFspTabOutputs.solutions.data.dps{j, ipar}.p,false);
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
%     legends = [legends char(app.NameTable.Data(iSpecies,2))];
end

%% Plot trajectories, each species have a color
speciesColors = {'b', 'r', 'g'};
figure()
for iSim = 1:3
    % Plot the trajectories
    for iSim = speciesToPlot
        stairs([0:length(sensmdist{iSim})], [sensmdist{iSim};0], speciesColors{iSim},'Linewidth',2);
        hold on
    end
end
%% Add Legend Title and Axes
legend(legends,'Location','best');
timeTitle = sprintf('Sensitivity at time: t = %1.2f ',T_array(j));
paramTitle = app.SensParDropDown.Value;
titleTitle = append(timeTitle,'w.r.t. ',paramTitle);
title(titleTitle);
xlabel('Species Count');
ylabel('Sensitivity');
set(gca,'Fontsize',20);
end


