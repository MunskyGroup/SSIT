function ModelMergeSortPars(app,SSITapp)
allPars = {};

if length(SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName)>1
    app.MergeModels.Enable=1;
else
    app.MergeModels.Enable=0;    
end

for i = 1:length(SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName)
    ModFile = [SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName{i},...
        SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName{i}];
    
    try
        load(ModFile,'Loc_app_str');
        allPars(end+1:end+length(Loc_app_str.ModelParameterTable.Data(:,1)),1) = Loc_app_str.ModelParameterTable.Data(:,1);
    catch
        load(ModFile,'Project');
        Loc_app_str = Project;
        allPars(end+1:end+length(Loc_app_str.ModelParameterTable.Data(:,1)),1) = Loc_app_str.ModelParameterTable.Data(:,1);
    end       
end

[allParsUnique] = unique(allPars);

K = zeros(1,length(allParsUnique),'logical');
for i=1:length(allParsUnique)
    if sum(strcmp(allPars,allParsUnique{i}))>1
        K(i)=1;
    end
end
allParsShared = allParsUnique(K);

app.MergeModelConstraints.Data = {'tmp',0};
app.MergeModelConstraints.Data(1:length(allParsShared),1) = allParsShared;
app.MergeModelConstraints.Data(:,2) = {inf};

app.ModelList.Value = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName;
app.ModelList.Value(end+1) = {'Models are not yet merged'};
end