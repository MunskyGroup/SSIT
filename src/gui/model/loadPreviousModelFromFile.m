function loadPreviousModelFromFile(app)
% Loads the model, the intial conditions, time span, and FSP constraints 
% into the GUI.
[array_file,array_path] = uigetfile('*.mat','Select Model')
if isequal(array_file,0)
    return
else
    addpath(array_path);
    
    Load_Dans_Model(app, array_file)
    updateModel(app)
    app.FspInitCondField.Value = '[2,0,0]';
    app.FspPrintTimesField.Value = '[0:0.5:14]';
    updateTimeSliderFsp(app);
    
    app.FspConstraintTable.Data = {'-x1';'-x2';'-x3';'x1';'x2';'x3';'x1+x2'};
    app.FspConstraintTable.Data(:,3) = {0;0;0;2;2;1400;2};
    app.FspConstraintTable.Data(:,2) = {'<';'<';'<';'<';'<';'<';'<'};
    
    app.tab_pars_Dan.Data = app.ModelParameterTable.Data;
    app.tab_inp_Dan.Data = app.ModelInputTable.Data;
    
    
    app.UITable6.Data(:,1) = 0:5:100;
    app.UITable6.Data(:,2) = 100*ones(21,1);
end
end
