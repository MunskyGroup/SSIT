function [] = updatePlotsInSensFsp(app)
% This function Updates the plots when any change occurs in selection of
% the species shown and changes in the timeslider.

% Find the time index to plot
T_array = unique(eval(app.SensPrintTimesEditField.Value));
[~,j] = min(abs(T_array-app.SensPlotTimeSlider.Value));

% Find the parameter index to plot
ipar = 1;
num_pars = length(app.ReactionsTabOutputs.parameters);
while (~strcmp(app.ReactionsTabOutputs.parameters{ipar}, app.SensParDropDown.Value) && ipar < num_pars)
    ipar = ipar + 1;
end

% Compute the marginal distributions
solutionFormat = app.SensFspTabOutputs.solutions.format;

Nd = app.SensFspTabOutputs.solutions.data{j}.p.dim;
for i=1:Nd
    INDS = setdiff([1:Nd],i);
    if ~isempty(INDS)
        mdist{i} = double(app.SensFspTabOutputs.solutions.data{j}.p.sumOver(INDS).data);
        sensmdist{i} = double(app.SensFspTabOutputs.solutions.data{j}.S(ipar).sumOver(INDS).data);
    else
        mdist{i} = double(app.SensFspTabOutputs.solutions.data{j}.p.data);
        sensmdist{i} = double(app.SensFspTabOutputs.solutions.data{j}.S(ipar).data);
    end
end

hold(app.SensProbAxes,'off'); %%% Add in bin-size allowance
hold(app.SensDerivativeAxes,'off'); %%% Add in bin-size allowance
N_max_x = 1e-6;
N_max_y = 1e-6;
% Plot the marginal distribution

for iSpecies = 1:Nd
    iPlot(iSpecies) = max(contains(app.SpeciesForSensPlot.Value,app.ReactionsTabOutputs.varNames{iSpecies}));
    if iPlot(iSpecies)
        stairs(app.SensProbAxes, [0:length(mdist{iSpecies})], [mdist{iSpecies};0],'Linewidth',2);
        hold(app.SensProbAxes,'on');
        N_max_x = max(N_max_x,length(mdist{iSpecies})-1);
        N_max_y = max(N_max_y,max(abs(mdist{iSpecies})));
    end
end

app.FspSolutionAtTimeLabel.Text = ['Solutions at time: t = ',num2str(T_array(j))];
app.SensProbAxes.XLim = [0,N_max_x + 1];
app.SensProbAxes.YLim = N_max_y*[0,1];
app.FSPXminEditField.Value = app.SensProbAxes.XLim(1);
app.FSPXmaxEditField.Value = app.SensProbAxes.XLim(2);
app.FSPYminEditField.Value = app.SensProbAxes.YLim(1);
app.FSPYmaxEditField.Value = app.SensProbAxes.YLim(2);

legend(app.SensProbAxes,app.ReactionsTabOutputs.varNames(iPlot),'Location','best')

% Plot the marginal sensitivity
N_max_y = 1e-6;
for iSpecies = 1:Nd
    if iPlot(iSpecies)
        stairs(app.SensDerivativeAxes, [0:length(sensmdist{iSpecies})-1], [sensmdist{iSpecies}],'Linewidth',2);
        hold(app.SensDerivativeAxes,'on');
        N_max_x = max(N_max_x,length(sensmdist{iSpecies})-1);
        N_max_y = max(N_max_y,max(abs(sensmdist{iSpecies})));
    end
end
app.FspSolutionAtTimeLabel.Text = ['Solutions at time: t = ',num2str(T_array(j))];
app.SensDerivativeAxes.XLim = [0,N_max_x + 1];
app.SensDerivativeAxes.YLim = N_max_y*[-1 1];

legend(app.SensDerivativeAxes,app.ReactionsTabOutputs.varNames(iPlot),'Location','best')
end