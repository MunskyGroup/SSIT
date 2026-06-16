function [] = updateModelswithinDropDown(app)
% Updates the models included under the selected model classification
% within the examples dropdown.

value = app.ModelDropDown.Value;
path(path,['Models/',app.ModelUsePresetExampleTypeDropDown.Value])
fileName = append('Models/',app.ModelUsePresetExampleTypeDropDown.Value,'/',value);
app.ModelFile.fileName = fileName;
J = strfind(value,'.');
app.ModelFile.modelName = value(1:J(end)-1);
if strcmp(value(end-2:end),'mat')
    [app] = loadModelBP(app, [], fileName);

    if ~exist(['RemovedModels/',app.ModelUsePresetExampleTypeDropDown.Value],'dir')
        mkdir(['RemovedModels/',app.ModelUsePresetExampleTypeDropDown.Value]);
    end
    movefile(fileName,['RemovedModels/',app.ModelUsePresetExampleTypeDropDown.Value]);
    
elseif strcmp(value(end-1:end),'.m')
    app.ReactionsTabOutputs.parameters={};
    app.ReactionsTabOutputs.presetParameters = {};
    ldmod = str2func(['@(x)',value(1:end-2),'(x)']);
    
    app.ReactionsTabOutputs.inputs={};
    app.ReactionsTabOutputs.presetInputs={};
    ldmod(app);
    
    app.ModelInputTable.Data = [app.ReactionsTabOutputs.inputs,app.ReactionsTabOutputs.presetInputs];       % Fills in the Input Table
    
    try
        app.ModelAbout.Value = app.ReactionsTabOutputs.modelInfo;
    catch
        app.ModelAbout.Value = {'About the Model';'';'Not provided'};
    end
    
    % Call code to update model and save .mat version for later use.
    updateModel(app,true,fileName,1);

    if ~exist(['RemovedModels/',app.ModelUsePresetExampleTypeDropDown.Value],'dir')
        mkdir(['RemovedModels/',app.ModelUsePresetExampleTypeDropDown.Value]);
    end
    movefile(fileName,['RemovedModels/',app.ModelUsePresetExampleTypeDropDown.Value]);
end
updateTimeSliderFsp(app);