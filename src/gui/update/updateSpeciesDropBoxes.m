% This script updates the species options for plotting to include distorted
% or non-distorted data.
speciesBoxes= {'SpeciestoShowListBoxMargFSP',...
    'SpeciestoShowListBox_2','SpeciestoShowListBoxMargFSPvT',...
    'SpeciestoShowListBoxMeans','JointSp1','JointSp2',...
    'SpeciesForSensPlot'};

species = app.SSITModel.species;
speciesStochastic = setdiff(species,app.SSITModel.hybridOptions.upstreamODEs);

if ~isempty(app.SSITModel.pdoOptions)&&isfield(app.SSITModel.pdoOptions,'PDO')&&~isempty(app.SSITModel.pdoOptions.PDO)
    if ~isfield(app.SSITModel.pdoOptions,'unobservedSpecies')
        observedSpecies = speciesStochastic;
        app.SSITModel.pdoOptions.unobservedSpecies = [];
    else
        observedSpecies = setdiff(speciesStochastic,app.SSITModel.pdoOptions.unobservedSpecies);
    end
    for i = 1:length(speciesBoxes)
        app.(speciesBoxes{i}).Items = ...
        [speciesStochastic;cellfun(@(s)[s ' (distorted)'], observedSpecies, 'UniformOutput', false)];
    end
    app.SpeciesDropDown.Items = setdiff(speciesStochastic,app.SSITModel.pdoOptions.unobservedSpecies);
    app.SpeciesDropDown.Visible = true;
    app.SolutionTimeDropDown.Visible = true;
    app.SetDistortionParametersButton.Enable = true;
    app.ShowDistortionPlotButton.Enable = true;
    app.SolutionTimeDropDown.Items = arrayfun(@num2str, app.SSITModel.tSpan, 'UniformOutput', 0);
    app.SolutionTimeDropDown.Value=app.SolutionTimeDropDown.Items(end);
else
    for i = 1:length(speciesBoxes)
        app.(speciesBoxes{i}).Items = speciesStochastic';
    end
    app.SpeciesDropDown.Visible = false;
    app.SolutionTimeDropDown.Visible = false;
    app.SetDistortionParametersButton.Enable = false;
    app.ShowDistortionPlotButton.Enable = false;
    app.PDO_Axis.Visible = false;
    app.PDO_Axis2.Visible = false;
end

%% Species boxes that do not remove the upstream ODE species when using hybrid models.
speciesBoxes = {'SpeciestoShowListBox'};
species = app.SSITModel.species;

if ~isempty(app.SSITModel.pdoOptions)&&isfield(app.SSITModel.pdoOptions,'PDO')&&~isempty(app.SSITModel.pdoOptions.PDO)
    speciesDistorted = setdiff(setdiff(species,app.SSITModel.pdoOptions.unobservedSpecies),app.SSITModel.hybridOptions.upstreamODEs);
else
    speciesDistorted = {};
end
%     if ~isfield(app.SSITModel.pdoOptions,'unobservedSpecies')
%         observedSpecies = speciesStochastic';
%         app.SSITModel.pdoOptions.unobservedSpecies = [];
%     else
%         observedSpecies = setdiff(speciesStochastic,app.SSITModel.pdoOptions.unobservedSpecies);
%     end
%     for i = 1:length(speciesBoxes)
%         app.(speciesBoxes{i}).Items = ...
%         [speciesStochastic;cellfun(@(s)[s ' (distorted)'], observedSpecies, 'UniformOutput', false)];
%     end
% else
for i = 1:length(speciesBoxes)
    app.(speciesBoxes{i}).Items = ...
        [species;cellfun(@(s)[s ' (distorted)'], speciesDistorted, 'UniformOutput', false)];
end
% end

%% Species boxes that will NOT include distorted data.
speciesBoxes = {'UpstreamSpeciesListBox'};
for i = 1:length(speciesBoxes)
    app.(speciesBoxes{i}).Items = app.SSITModel.species';
end

% Species boxes that should have no selections
if app.SSITModel.useHybrid
    app.UseHybridModelSwitch.Value = 'On';
    app.UpstreamSpeciesListBox.Visible = true;
    if ~isempty(app.SSITModel.hybridOptions.upstreamODEs)
        app.UpstreamSpeciesListBox.Value = app.SSITModel.hybridOptions.upstreamODEs;
    else
        app.UpstreamSpeciesListBox.Value = {};
    end
end

%% Species boxes that will NOT include upstream reactions.
speciesBoxes = {'ObservableSpeciesListBox',...
    'ObservableSpeciesListBox_2'};
for i = 1:length(speciesBoxes)
    app.(speciesBoxes{i}).Items = speciesStochastic;
    app.(speciesBoxes{i}).Value = setdiff(speciesStochastic,app.SSITModel.hybridOptions.upstreamODEs);
end

%% Species boxes that should have no selections
if app.SSITModel.useHybrid
    app.UseHybridModelSwitch.Value = 'On';
    app.UpstreamSpeciesListBox.Visible = true;
    if ~isempty(app.SSITModel.hybridOptions.upstreamODEs)
        app.UpstreamSpeciesListBox.Value = app.SSITModel.hybridOptions.upstreamODEs;
    else
        app.UpstreamSpeciesListBox.Value = {};
    end
else
    app.UseHybridModelSwitch.Value = 'Off';
    app.UpstreamSpeciesListBox.Visible = false;
    app.SSITModel.hybridOptions.upstreamODEs = {};
    app.UpstreamSpeciesListBox.Items = app.SSITModel.species;
    app.UpstreamSpeciesListBox.Value = {};
end


