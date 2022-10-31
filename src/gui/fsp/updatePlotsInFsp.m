function [] = updatePlotsInFsp(app)
% This function Updates the plots when any change occurs in selection of
% the species shown and changes in the timeslider.
T_array = eval(app.FspPrintTimesField.Value);
[~,j] = min(abs(T_array-app.FspTimeSlider.Value));

% Compute the marginal distributions
% mdist = ssit.fsp.marginals(app.FspTabOutputs.solutions{j}.states, app.FspTabOutputs.solutions{j}.p);
Nd = app.FspTabOutputs.solutions{j}.p.dim;
if Nd==1
    mdist{1} = double(app.FspTabOutputs.solutions{j}.p.data);
else
    for i=1:Nd
        INDS = setdiff([1:Nd],i);
        mdist{i} = double(app.FspTabOutputs.solutions{j}.p.sumOver(INDS).data);
    end
end

hold(app.FspAxes,'off'); %%% Add in bin-size allowance
N_max_x = 0;
N_max_y = 1e-6;

%% Check What species to plot
species2Plot=[];
legends={};
nSpecies = length(app.ReactionsTabOutputs.varNames);
for iSpecies = 1:nSpecies
    if max(contains(app.SpeciestoShowListBox_2.Value,app.ReactionsTabOutputs.varNames{iSpecies}))
        species2Plot = [species2Plot iSpecies];
        legends=[legends char(app.NameTable.Data(iSpecies,2))];
    end
end
% Plot FSP graph for each selected species
hold(app.FspAxes,'off');
speciesColors = {'b', 'r', 'g', 'k', 'm', 'c', 'b--', 'r--', 'g--', 'k--', 'm--', 'c--'};

% Plot the trajectories
for iSpecies = species2Plot
    stairs(app.FspAxes, [0:length(mdist{iSpecies})], [mdist{iSpecies};0],speciesColors{iSpecies} ,'Linewidth',2);
    hold(app.FspAxes,'on');
    N_max_x = max(N_max_x,length(mdist{iSpecies})-1);
    N_max_y = max(N_max_y,max(mdist{iSpecies}));
end

% Adds legend
legend(app.FspAxes,legends,'Location','northeast')

app.FspSolutionAtTimeLabel.Text = ['Solutions at time: t = ',num2str(T_array(j))];
app.FspAxes.XLim = [0,N_max_x + 1];
app.FspAxes.YLim = [0,N_max_y];
app.FSPXminEditField.Value = app.FspAxes.XLim(1);
app.FSPXmaxEditField.Value = app.FspAxes.XLim(2);
app.FSPYminEditField.Value = app.FspAxes.YLim(1);
app.FSPYmaxEditField.Value = app.FspAxes.YLim(2);

end
