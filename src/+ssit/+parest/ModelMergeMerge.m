function ModelMergeMerge(app,SSITapp,saveIt,fspBounds,overWriteParameters)
arguments
    app
    SSITapp
    saveIt = true;
    fspBounds = [];
    overWriteParameters = true;
end

nMods = length(SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName);

allParsRepeated = app.MergeModelConstraints.Data(:,1);
parSimilarityConstraints = [app.MergeModelConstraints.Data{:,2}];
allParsStrict = allParsRepeated(isinf(parSimilarityConstraints));
allParsConstr = allParsRepeated(~isinf(parSimilarityConstraints));
parSimilarityConstraints = parSimilarityConstraints(~isinf(parSimilarityConstraints));

% Load all models in in list and make a list of names and initial values.
ModelData = {};
for i = 1:nMods
    ModFile = [SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName{i},...
        SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName{i}];

    try
        LAS{i} = load(ModFile,'Loc_app_str','vars_to_fit');
        modPar{i}.Names = LAS{i}.Loc_app_str.ModelParameterTable.Data(:,1);
    catch
        load(ModFile,'Project');
        LAS{i}.Loc_app_str = Project;
        LAS{i}.vars_to_fit = find(strcmp(Project.fit_parameters_table.Data(:,3),'y'));
        modPar{i}.Names = LAS{i}.Loc_app_str.ModelParameterTable.Data(:,1);
    end

    % Load names of parameters in the current model
    modPar{i}.Values = [LAS{i}.Loc_app_str.ModelParameterTable.Data{:,2}];
    if app.SuppressFSPExpansion.Value
        LAS{i}.Loc_app_str.FspErrorTolField.Value = inf;
    end
    
    % Update FSP bounds from previous fit or FSP eval.
    if ~isempty(fspBounds)
        LAS{i}.Loc_app_str.FspTabOutputs.bounds = fspBounds{i};
        for j = 1:size(LAS{i}.Loc_app_str.FspConstraintTable.Data,1)
            LAS{i}.Loc_app_str.FspConstraintTable.Data{j,3} = fspBounds{i}(j);
        end
    end
    ModelData{i} = LAS{i}.Loc_app_str;
end
SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.ModelData = ModelData;

% Generate mapping between lists parameters for each model (Local) and non-redudant
% list of parameters for all models (Global) and a linear Map to generate
% the local parameters from the global parameter vector.
kGlob = length(allParsStrict);
MapAllToEach = [];
StrictParsCell = cell(length(allParsStrict),1);
MapConstraints = cell(length(allParsConstr),1);
clear OverallParGuesses OverallParNames

for i=1:nMods
    kLoc = 0;
    for j=1:length(modPar{i}.Names)
        kLoc = kLoc+1;
        switch modPar{i}.Names{j}
            case allParsStrict % strict constraint use lumped model.
                kk = find(strcmp(allParsStrict,modPar{i}.Names{j}));
                StrictParsCell(kk) = {[StrictParsCell{kk},modPar{i}.Values(j)]};
            case allParsConstr % constraint
                kc = find(strcmp(allParsConstr,modPar{i}.Names{j}));
                kGlob = kGlob+1; kk=kGlob;
                MapConstraints(kc) = {[MapConstraints{kc},kGlob]};
                OverallParGuesses(kGlob,1) = modPar{i}.Values(j);
                OverallParNames(kGlob) = {[modPar{i}.Names{j},'_M',num2str(i)]};
            otherwise % completely free
                kGlob = kGlob+1; kk=kGlob;
                OverallParGuesses(kGlob,1) = modPar{i}.Values(j);
                OverallParNames(kGlob) = {[modPar{i}.Names{j},'_M',num2str(i)]};
        end
        MapAllToEach{i}(kLoc,kk) = 1;       
    end
end

% Fix sizes of Maps to match total number of global parameters
for i = 1:nMods
    if size(MapAllToEach{i},2)<kGlob
        MapAllToEach{i}(1,kGlob)=0;
    end
end
SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.MapAllToEach = MapAllToEach;

% Replace strict constrained parameters with mean of their guesses.
for i=1:length(allParsStrict)
    OverallParGuesses(i,1) = mean(StrictParsCell{i});
    OverallParNames(i) = allParsStrict(i);
end

% Update fit parameter table with the new parameter names.
if overWriteParameters
    app.fit_parameters_table.Data = {};
    for i=1:length(OverallParNames)
        app.fit_parameters_table.Data(i,1:3) = {OverallParNames{i},OverallParGuesses(i,1),'y'};
    end
    
    % Adjust initial guess parameters to log space for ease in fitting.
    OverallParGuesses = log10(OverallParGuesses);
end

% Formulate objective functions for each model/data combination.
for i=1:nMods
    OBJ{i} =  @(x)-get_fit_error(10.^[MapAllToEach{i}*x],LAS{i}.Loc_app_str,LAS{i}.vars_to_fit);
end

% Define function for constrain penalty calculation.
OBJ_Constraint = @(x)ComputeConstraintViolation(10.^x,MapConstraints,parSimilarityConstraints);

% Combine all the objective functions. 
OverallObjective = @(x) ( getOverallObj(OBJ,x) + OBJ_Constraint(x) );

% Define function to allow for plotting results.
for i=1:nMods
    FitResults{i} = @(x) returnFitResults(10.^[MapAllToEach{i}*x] , LAS{i}.Loc_app_str , LAS{i}.vars_to_fit );
end


SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FitResults = FitResults;
SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.OverallObjective = OverallObjective;
SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.OverallParGuesses = OverallParGuesses;

app.MergedModelDropDown.Enable = 1;
app.MergedModelDropDown.Items = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName;
app.MergedModelDropDown.Value = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName{1};

app.MergedModelDropDown2.Enable = 1;
app.MergedModelDropDown2.Items = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName;
app.MergedModelDropDown2.Value = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName{end};

if saveIt
    ssit.parest.SaveProject(app,SSITapp)
    figure(app.UIFigure);
end

ssit.parest.makeMergeModelPlots(app,app.SSITapp,true); 

end

function Total = getOverallObj(OBJ,x)
Total = 0;
for i = 1:length(OBJ)
    Total = Total + OBJ{i}(x);
end
end

function [constrPenalty] = ComputeConstraintViolation(x,MapConstraints,parSimilarityConstraints)
constrPenalty = 0;
for i = 1:length(parSimilarityConstraints)
    if isfinite(parSimilarityConstraints(i))&&parSimilarityConstraints(i)~=0
        constrPenalty = constrPenalty + ...
            parSimilarityConstraints(i)*sum((x(MapConstraints{i})-mean(x(MapConstraints{i}))).^2);
    end
end
end

function [fit_error,app] = get_fit_error(x,app,vars_to_fit)
app.fit_parameters_table.Data(vars_to_fit,2) = num2cell(x);
app.ReactionsTabOutputs.parameters(:,2) = app.fit_parameters_table.Data(:,2);
app.ModelParameterTable.Data=app.fit_parameters_table.Data(:,1:2);
app = ssit.parest.updateModelSolveAndCompareToData(app);
fit_error = app.DataLoadingAndFittingTabOutputs.J_LogLk;
end

function [fitResults] = returnFitResults(x,app,vars_to_fit)
app.fit_parameters_table.Data(vars_to_fit,2) = num2cell(x);
app.ModelParameterTable.Data=app.fit_parameters_table.Data(:,1:2);
app.ReactionsTabOutputs.parameters(:,2) = app.fit_parameters_table.Data(:,2);
app = ssit.parest.updateModelSolveAndCompareToData(app);
fitResults = {app.DataLoadingAndFittingTabOutputs.fitResults.current,...
    app.DataLoadingAndFittingTabOutputs.fitResults.currentData,...
    app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,...
    app.FspTabOutputs.bounds};
end


