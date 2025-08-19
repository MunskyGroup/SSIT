function app = ChooseSSITModelValue(app)
Model = load(app.ModelFile.fileName,app.ChooseSSITModel.Value);
app.SSITModel = Model.(app.ChooseSSITModel.Value);
app.ModelFile.modelName = app.ChooseSSITModel.Value;
app.ModelFile.propFileName = []; % Clear name of propensity functions.
updateAppFromSSIT(app)

%% Update other GUI Fields from SSIT Model
if ~isempty(app.SSITModel.GUIProps)
    for i = 1:length(app.SSITModel.GUIProps.SaveItems)
        app.(app.SSITModel.GUIProps.SaveItems{i}).Items = app.SSITModel.GUIProps.Items.(app.SSITModel.GUIProps.SaveItems{i});
    end

    for i = 1:length(app.SSITModel.GUIProps.SaveValue)
        app.(app.SSITModel.GUIProps.SaveValue{i}).Value = app.SSITModel.GUIProps.Values.(app.SSITModel.GUIProps.SaveValue{i});
    end

    % for i = 1:length(app.SSITModel.GUIProps.SaveText)
    %     app.(app.SSITModel.GUIProps.SaveValue{i}).Text = app.SSITModel.GUIProps.Texts.(app.SSITModel.GUIProps.SaveText{i});
    % end

    for i = 1:length(app.SSITModel.GUIProps.SaveData)
        app.(app.SSITModel.GUIProps.SaveData{i}).Data = app.SSITModel.GUIProps.Data.(app.SSITModel.GUIProps.SaveData{i});
    end

    for i = 1:length(app.SSITModel.GUIProps.SaveFields)
        app.(app.SSITModel.GUIProps.SaveFields{i}) = app.SSITModel.GUIProps.Fields.(app.SSITModel.GUIProps.SaveFields{i});
    end

else
    clearAppData(app);
end

info = dir(app.ModelFile.fileName);

app.FileModelLabel.Text = {['File: ',app.ModelFile.fileName];...
    ['Model: ',app.ModelFile.modelName];...
    ['Last Saved: ',info.date]};