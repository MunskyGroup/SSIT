function [] = plotSsaInGui(app)
% This function is meant to plot the SSA information into the GUI for the
% user to evaluate

%% Extract SSA samples from storage
if isempty(app.SSITModel.Solutions)||~isfield(app.SSITModel.Solutions,'trajs')
    msgbox('There is no trajectory to plot. Please click runSsa first to generate trajectories.');
    return
end

samples = app.SSITModel.Solutions.trajs;
if isfield(app.SSITModel.Solutions,'trajsDistorted')
    speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);
    for iSp = 1:length(speciesStochastic)
        if ~max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic{iSp}))
            samples(end+1,:,:) = app.SSITModel.Solutions.trajsDistorted(iSp,:,:);
        end
    end
end

%% Extract the number and timepoints to plot
if app.SsaNumSimField.Value<app.SsaNumSimToPlotField.Value
    app.SsaNumSimToPlotField.Value = app.SsaNumSimField.Value;
end
numToPlot = app.SsaNumSimToPlotField.Value;  % Pulls the number of simulations
numTotal = app.SsaNumSimField.Value;  % Pulls the total number of samples available

if isempty(samples)
    msgbox('There is no trajectory to plot. Please click runSsa first to generate trajectories.');
end
timesToPlot = eval(app.PrintTimesEditField.Value);  % Pulls the time array
%% Check what species to plot
speciesToPlot = [];
legends = {};
nSpecies = length(app.SpeciestoShowListBox.Items);
for iSpecies = 1:nSpecies
    if max(strcmpi(app.SpeciestoShowListBox.Value,app.SpeciestoShowListBox.Items{iSpecies}))
        speciesToPlot = [speciesToPlot iSpecies];
        legends = [legends app.SpeciestoShowListBox.Items(iSpecies)];
    end
end
%% Plot trajectories, each species have a color
LG = legends;
lg = {};
ode = [];
hold(app.SsaTrajAxes,'off');
speciesColors = {'b', 'r', 'g', 'm', 'c', 'b--', 'r--', 'g--', 'm--', 'c--'};
speciesColors2 = {'k', 'k--', 'k.-', 'k-x', 'k-s', 'k-o'};
for isim = 1:numToPlot
    % Plot the trajectories
    for iSpecies = speciesToPlot
        trajectories(iSpecies,:,isim)=plot(app.SsaTrajAxes,timesToPlot, samples(iSpecies,:,isim), ...
            speciesColors{iSpecies}); hold(app.SsaTrajAxes,'on');
    end
end

for iSpecies = speciesToPlot(speciesToPlot<=length(app.SSITModel.species))
    if app.SsaShowOdeCheckBox.Value == 1
        if isempty(app.SSITModel.Solutions.ode)
            runOde(app)
        end
        ode(iSpecies) = plot(app.SsaTrajAxes,app.SSITModel.Solutions.T_array,app.SSITModel.Solutions.ode(:,iSpecies),speciesColors2{iSpecies},'LineWidth',3); hold(app.SsaTrajAxes,'on');
        LG{end+1} = [append(char(app.SSITModel.species{iSpecies}),' ode')];
    end
end
% Find plottable graphics
odeInd = find(ode);
for iOde = 1:length(odeInd)
    odePlot(iOde) = ode(odeInd(iOde));
end
if isempty(ode)
    graphicArray = [trajectories(speciesToPlot(1):speciesToPlot(end))];
else
    graphicArray = [trajectories(speciesToPlot(1):speciesToPlot(end)) odePlot];
end
arrayTest = isgraphics(graphicArray);
arrayInd = find(arrayTest); %find graphical indicies for legend
graphicArrayPlot = gobjects(1,length(arrayInd)); %Initialize graphical object
for iG = 1:length(arrayInd)
    graphicArrayPlot(iG) = graphicArray(arrayInd(iG));
end
% Add legends
legend(app.SsaTrajAxes,[graphicArrayPlot],LG , 'Location', 'best');

%% Histogram plot
[~, timeIndex] = min(abs(timesToPlot-app.SsaTimeSlider.Value));
maxN = max(samples(speciesToPlot,timeIndex,:),[],'all');
minN = min(samples(speciesToPlot,timeIndex,:),[],'all');
binSize = max(1,ceil((maxN-minN)/sqrt(numTotal)));
bins = [minN:binSize:maxN+binSize];
if numTotal>1
    hold(app.SsaHistAxes,'off');
    for iSpecies = speciesToPlot
        [h] = histcounts(squeeze(samples(iSpecies,timeIndex,:)),bins,'Normalization','probability');
        stairs(app.SsaHistAxes,bins(1:end-1),h,speciesColors{iSpecies},'LineWidth',2);  hold(app.SsaHistAxes,'on');
    end
    title(app.SsaHistAxes,['Histogram at t = ',num2str(timesToPlot(timeIndex))]);
else
    cla(app.SsaHistAxes); % Clears the axes when there are no simulations available
end
app.SsaSolutionsattimeLabel.Text = ['Solutions at time: t = ',num2str(timesToPlot(timeIndex))];
app.SSA_hist_XminEditField.Value = app.SsaHistAxes.XLim(1);
app.SSA_hist_XmaxEditField.Value = app.SsaHistAxes.XLim(2);
app.SSA_hist_YminEditField.Value = app.SsaHistAxes.YLim(1);
app.SSA_hist_YmaxEditField.Value = app.SsaHistAxes.YLim(2);
% Add legends
legend(app.SsaHistAxes, legends, 'Location', 'best');
end
