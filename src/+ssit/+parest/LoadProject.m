function LoadProject(app,SSITapp)
[fn,pn] = uigetfile('*.mat','Open Previous Project:');
load([pn,fn],'ProjectA','ProjectB','SSITFields','AppFields')
app.ProjectFileInfo.pn = pn;
app.ProjectFileInfo.fn = fn;
app.ProjectFileInfo.ProjName = fn(1:end-4);

for i=1:length(SSITFields)
    try
    eval(['SSITapp.',SSITFields{i},' = ProjectA.',SSITFields{i},';'])
    catch
    end
end
for i=1:length(AppFields)
    try
        eval(['app.',AppFields{i},' = ProjectB.',AppFields{i},';'])
    catch
    end
end

app.ModelList.Value = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName;
app.ModelList.Value(end+1:end+2) = {'Merged Model Saved as:',fn};
app.MergeModels.Enable = 1;
app.ProjectName.Value = fn(1:end-4);

app.MergedModelDropDown.Enable = 1;
app.MergedModelDropDown.Items = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName;
app.MergedModelDropDown.Value = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName{1};

app.MergedModelDropDown2.Enable = 1;
app.MergedModelDropDown2.Items = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName;
app.MergedModelDropDown2.Value = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName{end};

figure(app.UIFigure);
