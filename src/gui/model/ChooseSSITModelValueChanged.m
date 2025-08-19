function app = ChooseSSITModelValueChanged(app)
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
        app.SSITModel.GUIProps.Values.(app.SSITModel.GUIProps.SaveValue{i}) = app.(app.SSITModel.GUIProps.SaveValue{i}).Value;
    end

    for i = 1:length(app.SSITModel.GUIProps.SaveData)
        app.(app.SSITModel.GUIProps.SaveData{i}).Data = app.SSITModel.GUIProps.Data.(app.SSITModel.GUIProps.SaveData{i});
    end

    for i = 1:length(app.SSITModel.GUIProps.SaveFields)
        app.(app.SSITModel.GUIProps.SaveFields{i}) = app.SSITModel.GUIProps.Fields.(app.SSITModel.GUIProps.SaveFields{i});
    end

    if ~isempty()
        app.SSITModel
    end

end