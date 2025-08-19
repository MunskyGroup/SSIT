function app =clearAppData(app)
% This function resets the SSITGUI app to remove choices when the model is
% changed.

app.SSITModel.dataSet = [];
app.FieldsinDataTextArea.Value = 'None - need to load data first.';
app.TotalCellsInDataLabel.Text = 'Total Cells In Data:';
app.NumberafterConstraintsLabel.Text = 'Number after Constraints:';

ListItems = {'ParEstFitTimesList','ObservableSpeciesListBox',...
    'DataSpecies1','DataSpecies2','DataSpecies3','DataSpecies4','DataSpecies5',...
    'DataSpecies6','DataSpecies7','DataSpecies8','DataSpecies9',...
    };

for i = 1:length(ListItems)
    app.(ListItems{i}).Items = {};
end

app = clearDataConstraints(app);

app.ObservableSpeciesListBox.Items = app.SSITModel.species;