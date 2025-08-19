function saveGUI(app)

if app.IncludeGUIDatainSaveSwitch.Value
    SaveItems = {'ParEstFitTimesList','ObservableSpeciesListBox',...
        'DataSpecies1','DataSpecies2','DataSpecies3','DataSpecies4','DataSpecies5','DataSpecies6','DataSpecies7','DataSpecies8','DataSpecies9',...
        'DataConstrChoice1','DataConstrChoice2','DataConstrChoice3','DataConstrChoice4','DataConstrChoice5','DataConstrChoice6',...
        'AndOr2','AndOr3','AndOr4','AndOr5','AndOr6',...
        'DataLogical1','DataLogical3','DataLogical4','DataLogical5','DataLogical6'};
    SaveValue = {'SsaNumSimField','PrintTimesEditField','PrintTimesEditField',...
        'FspPrintTimesField','FspPiecewiseCheckBox','FspAutoExpandCheckBox',...
        'DataSpecies1','DataSpecies2','DataSpecies3','DataSpecies4','DataSpecies5','DataSpecies6','DataSpecies7','DataSpecies8','DataSpecies9',...
        'DataConstrChoice1','DataConstrChoice2','DataConstrChoice3','DataConstrChoice4','DataConstrChoice5','DataConstrChoice6',...
        'AndOr2','AndOr3','AndOr4','AndOr5','AndOr6',...
        'DataLogical1','DataLogical3','DataLogical4','DataLogical5','DataLogical6',...
        'DataConstrText1','DataConstrText2','DataConstrText3','DataConstrText4','DataConstrText5','DataConstrText6',...
        'ParEstFitTimesList','ObservableSpeciesListBox','FieldsinDataTextArea'};
    SaveText = {'Species1Label','Species2Label','Species3Label','Species4Label','Species5Label','Species6Label','Species7Label','Species8Label','Species9Label',...
        'TotalCellsInDataLabel','NumberafterConstraintsLabel'};
    SaveData = {'FspConstraintTable','SensParameterSelectionTable',...
        'fit_parameters_table'};
    SaveFields = {'DataLoadingAndFittingTabOutputs'};

    app.SSITModel.GUIProps.SaveItems = SaveItems;
    app.SSITModel.GUIProps.SaveValue = SaveValue;
    app.SSITModel.GUIProps.SaveText = SaveText;
    app.SSITModel.GUIProps.SaveData = SaveData;
    app.SSITModel.GUIProps.SaveFields = SaveFields;
    
    for i = 1:length(SaveItems)
        app.SSITModel.GUIProps.Items.(SaveItems{i}) = app.(SaveItems{i}).Items;
    end

    for i = 1:length(SaveValue)
        app.SSITModel.GUIProps.Values.(SaveValue{i}) = app.(SaveValue{i}).Value;
    end

    for i = 1:length(SaveText)
        app.SSITModel.GUIProps.Texts.(SaveText{i}) = app.(SaveText{i}).Text;
    end

    for i = 1:length(SaveData)
        app.SSITModel.GUIProps.Data.(SaveData{i}) = app.(SaveData{i}).Data;
    end

    for i = 1:length(SaveFields)
        app.SSITModel.GUIProps.Fields.(SaveFields{i}) = app.(SaveFields{i});
    end

    
end

