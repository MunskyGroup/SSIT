function app =clearAppData(app)
% This function resets the SSITGUI app to remove choices when the model is
% changed.

if ~isempty(app.SSITModel.dataSet)
    app.FieldsinDataTextArea.Value = 'Not available for loaded model.';
    app.TotalCellsInDataLabel.Text = ['Total Cells In Data:: ',num2str(size(app.SSITModel.dataSet.DATA,1))];
    app.NumberafterConstraintsLabel.Text = ['Number after Constraints: ',num2str(sum(app.SSITModel.dataSet.nCells))];
    app.DataFileNoneSelectedLabel.Text = 'Data File: Unknown';
    
else
    app.FieldsinDataTextArea.Value = 'None - need to load data first.';
    app.TotalCellsInDataLabel.Text = 'Total Cells In Data:';
    app.NumberafterConstraintsLabel.Text = 'Number after Constraints:';
    app.DataFileNoneSelectedLabel.Text = 'Data File: <None Selected>';
end

ListItems = {'ParEstFitTimesList','ObservableSpeciesListBox',...
    'DataSpecies1','DataSpecies2','DataSpecies3','DataSpecies4','DataSpecies5',...
    'DataSpecies6','DataSpecies7','DataSpecies8','DataSpecies9',...
    };

for i = 1:length(ListItems)
    app.(ListItems{i}).Items = {};
end

if ~isempty(app.SSITModel.dataSet)
    app.ParEstFitTimesList.Items = cellfun(@(s)[num2str(s)], num2cell(app.SSITModel.dataSet.times), 'UniformOutput', false);
    app.ParEstFitTimesList.Value = cellfun(@(s)[num2str(s)], num2cell(app.SSITModel.dataSet.times), 'UniformOutput', false);
end

app = clearDataConstraints(app);

app.ObservableSpeciesListBox.Items = app.SSITModel.species;