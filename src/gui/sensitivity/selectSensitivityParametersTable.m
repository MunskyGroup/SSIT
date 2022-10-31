function [] = selectSensitivityParametersTable(app)
% This function is meant to auto-populate the Sensitivity Parameters tab
% with edits in the reaction tab within the FSP GUI.

TMP(:,1) = app.ReactionsTabOutputs.parameters(:,1);
TMP(:,2) = app.ReactionsTabOutputs.parameters(:,2);

[j,~] = size(TMP); % gives the number of parameters

for i = 1:j
    TMP(i,3) = {'y'};
end

app.SensParameterSelectionTable.Data = TMP;
% app.SensInputTable.Data = app.ReactionsTabOutputs.inputs;

% app.tab_pars_Dan.Data = app.ModelParameterTable.Data;
% app.tab_inp_Dan.Data = app.ModelInputTable.Data;

app.fit_parameters_table.Data = [];
for i = 1:size(app.ModelParameterTable.Data,1)
    app.fit_parameters_table.Data{i,1} = app.ModelParameterTable.Data{i,1};
    app.fit_parameters_table.Data{i,2} = app.ModelParameterTable.Data{i,2};
    app.fit_parameters_table.Data{i,3} = 'y';
end
% app.fit_inputs_table.Data = app.ModelInputTable.Data;

end
