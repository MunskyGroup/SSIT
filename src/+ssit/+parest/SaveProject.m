function SaveProject(app,SSITapp,pn,fn)
arguments
    app
    SSITapp
    pn = [];
    fn = [];
end
SSITFields = {'ReactionsTabOutputs',...
    'DataLoadingAndFittingTabOutputs',...
    'ModelParameterTable.Data'...
    'ModelReactionTable.Data',...
    'NameTable.Data',...
    'ModelInputTable.Data'};

AppFields = {'MergeModelConstraints.Data',...
    'ModelList.Value',...
    'fit_parameters_table.Data',...
    'MergedModelDropDown.Items'};

for i=1:length(SSITFields)
    try
    eval(['ProjectA.',SSITFields{i},' = SSITapp.',SSITFields{i},';'])
    catch
    end
end
for i=1:length(AppFields)
    try
    eval(['ProjectB.',AppFields{i},' = app.',AppFields{i},';'])
    catch
    end
end

if isempty(fn)
    [fn,pn] = uiputfile('*.mat','Save Project as',app.ProjectFileInfo.ProjName);
    app.ProjectFileInfo.pn = pn;
    app.ProjectFileInfo.fn = fn;
    app.ProjectFileInfo.ProjName = fn(1:end-4);
end
save([pn,fn],'ProjectA','ProjectB','SSITFields','AppFields')

app.ModelList.Value = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName;
app.ModelList.Value(end+1:end+2) = {'Merged Model Saved as:',fn};
app.ProjectName.Value = fn(1:end-4);
