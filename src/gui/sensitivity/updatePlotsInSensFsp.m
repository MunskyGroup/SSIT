function [] = updatePlotsInSensFsp(app)
% This function Updates the plots when any change occurs in selection of
% the species shown and changes in the timeslider.

if isempty(app.SSITModel.Solutions)||~isfield(app.SSITModel.Solutions,'sens')
    app.SSITModel.solutionScheme = 'fspsens';
    [~,~,app.SSITModel] = app.SSITModel.solve;
end

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
Nd = app.SSITModel.Solutions.sens.data{j}.p.dim;
for i=1:Nd
    INDS = setdiff([1:Nd],i);
    if ~isempty(INDS)
        mdist{i} = double(app.SSITModel.Solutions.sens.data{j}.p.sumOver(INDS).data);
        sensmdist{i} = double(app.SSITModel.Solutions.sens.data{j}.S(ipar).sumOver(INDS).data);
    else
        mdist{i} = double(app.SSITModel.Solutions.sens.data{j}.p.data);
        sensmdist{i} = double(app.SSITModel.Solutions.sens.data{j}.S(ipar).data);
    end
end

%% Check What species to plot
species2Plot=[];
legends={};
nSpecies = length(app.SpeciesForSensPlot.Items);
for iSpecies = 1:nSpecies
    if max(strcmpi(app.SpeciesForSensPlot.Value,app.SpeciesForSensPlot.Items{iSpecies}))
        species2Plot = [species2Plot iSpecies];
        legends=[legends char(app.SpeciesForSensPlot.Items{iSpecies})];
    end
end

speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);

% Regenerate PDO if needed.
if isfield(app.SSITModel.pdoOptions,'PDO')&&max(species2Plot)>length(speciesStochastic)
    maxNum = app.SSITModel.Solutions.sens.data{j}.p.data.size;
    kSp = 0;
    for iS = 1:length(speciesStochastic)
        if max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic(iS)))
            maxNum(iS) = 0;
            curNum(iS) = 0;
        else
            kSp = kSp+1;
            curNum(iS) = size(app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp},2);
        end
    end
    if max(maxNum-curNum)>0
        [~,app.SSITModel] = app.SSITModel.generatePDO([],[],[],[],maxNum);
    end

    kSp = 0;
    for iSp = 1:length(speciesStochastic)
        if ~max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic{iSp}))
            INDS = setdiff([1:Nd],iSp);
            if ~isempty(INDS)
                px = double(app.SSITModel.Solutions.sens.data{j}.p.sumOver(INDS).data);
                Sx = double(app.SSITModel.Solutions.sens.data{j}.S(ipar).sumOver(INDS).data);
            else
                px = double(app.SSITModel.Solutions.sens.data{j}.p.data);
                Sx = double(app.SSITModel.Solutions.sens.data{j}.S(ipar).data);
            end
            kSp = kSp+1;
            mdist{end+1} = app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp}(:,1:length(px))*px;
            sensmdist{end+1} = app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp}(:,1:length(Sx))*Sx;
        end
    end
end

hold(app.SensProbAxes,'off'); %%% Add in bin-size allowance
hold(app.SensDerivativeAxes,'off'); %%% Add in bin-size allowance
N_max_x = 1e-6;
N_max_y = 1e-6;
% Plot the marginal distribution

for iSpecies = species2Plot
    stairs(app.SensProbAxes, [0:length(mdist{iSpecies})], [mdist{iSpecies};0],'Linewidth',2);
    hold(app.SensProbAxes,'on');
    N_max_x = max(N_max_x,length(mdist{iSpecies})-1);
    N_max_y = max(N_max_y,max(abs(mdist{iSpecies})));
end

% Adds legend
legend(app.SensProbAxes,legends,'Location','northeast')


app.FspSolutionAtTimeLabel.Text = ['Solutions at time: t = ',num2str(T_array(j))];
app.SensProbAxes.XLim = [0,N_max_x + 1];
app.SensProbAxes.YLim = N_max_y*[0,1];
app.FSPXminEditField.Value = app.SensProbAxes.XLim(1);
app.FSPXmaxEditField.Value = app.SensProbAxes.XLim(2);
app.FSPYminEditField.Value = app.SensProbAxes.YLim(1);
app.FSPYmaxEditField.Value = app.SensProbAxes.YLim(2);

% legend(app.SensProbAxes,app.SSITModel.species(iPlot),'Location','best')

% Plot the marginal sensitivity
N_max_y = 1e-6;
for iSpecies = species2Plot
    stairs(app.SensDerivativeAxes, [0:length(sensmdist{iSpecies})-1], [sensmdist{iSpecies}],'Linewidth',2);
    hold(app.SensDerivativeAxes,'on');
    N_max_x = max(N_max_x,length(sensmdist{iSpecies})-1);
    N_max_y = max(N_max_y,max(abs(sensmdist{iSpecies})));
end
app.FspSolutionAtTimeLabel.Text = ['Solutions at time: t = ',num2str(T_array(j))];
app.SensDerivativeAxes.XLim = [0,N_max_x + 1];
app.SensDerivativeAxes.YLim = N_max_y*[-1 1];

% Adds legend
legend(app.SensProbAxes,legends,'Location','northeast')
% legend(app.SensDerivativeAxes,app.SSITModel.species(iPlot),'Location','best')
end