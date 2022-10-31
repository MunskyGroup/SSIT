function ModelMergeAdd(app,SSITapp)
[FILENAME,PATHNAME] = uigetfile('*.mat','MultiSelect','on');
if iscell(FILENAME)
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName(end+1:end+length(FIELNAME)) = FILENAME;
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName(end+1:end+length(FIELNAME)) = {PATHNAME};
else
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName(end+1) = {FILENAME};
    SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName(end+1) = {PATHNAME};
end
ssit.parest.ModelMergeSortPars(app,SSITapp);
end
