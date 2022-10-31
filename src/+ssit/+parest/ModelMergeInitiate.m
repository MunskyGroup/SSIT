function ModelMergeInitiate(app,SSITapp)           
% Ask user to choose first model.
[FILENAME, PATHNAME] = uigetfile('*.mat','MultiSelect','on');  
SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName = {};
SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName = {};
if iscell(FILENAME)
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName = FILENAME;
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName(1:length(FILENAME)) = {PATHNAME};
else
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName = {FILENAME};
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName = {PATHNAME};
end

ssit.parest.ModelMergeSortPars(app,SSITapp);
figure(app.UIFigure);
end
