% This script updates the species options for plotting to include distorted
% or non-distorted data.
speciesBoxes = {'SpeciestoShowListBoxMargFSP','SpeciestoShowListBox',...
    'SpeciestoShowListBox_2','SpeciestoShowListBoxMargFSPvT',...
    'SpeciestoShowListBoxMeans','JointSp1','JointSp2',...
    'SpeciesForSensPlot',};
if ~isempty(app.SSITModel.pdoOptions)&&isfield(app.SSITModel.pdoOptions,'PDO')&&~isempty(app.SSITModel.pdoOptions.PDO)
    if ~isfield(app.SSITModel.pdoOptions,'unobservedSpecies')
        observedSpecies = app.SSITModel.species';
        app.SSITModel.pdoOptions.unobservedSpecies = [];
    else
        observedSpecies = setdiff(app.SSITModel.species,app.SSITModel.pdoOptions.unobservedSpecies);
    end
    for i = 1:length(speciesBoxes)
        app.(speciesBoxes{i}).Items = ...
        [app.SSITModel.species',cellfun(@(s)[s ' (distorted)'], observedSpecies, 'UniformOutput', false)];
    end
else
    for i = 1:length(speciesBoxes)
        app.(speciesBoxes{i}).Items = app.SSITModel.species';
    end
end
