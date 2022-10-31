function ModelMergeMerge2(app,SSITapp)

nMods = length(SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName);

allParsRepeated = app.MergeModelConstraints.Data(:,1);
parSimilarityConstraints = [app.MergeModelConstraints.Data{:,2}];
allParsStrict = allParsRepeated(isinf(parSimilarityConstraints));
allParsConstr = allParsRepeated(~isinf(parSimilarityConstraints));

% Load all models in in list and make a list of names and initial values.
for i = 1:nMods
    ModFile = [SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.PathName{i},...
        SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FileName{i}];
    LAS{i} = load(ModFile,'Loc_app_str','vars_to_fit');
    
    % Load names of parameters in the current model
    modPar{i}.Names = LAS{i}.Loc_app_str.ModelParameterTable.Data(:,1);
    modPar{i}.Values = [LAS{i}.Loc_app_str.ModelParameterTable.Data{:,2}];
end

% Generate mapping between lists parameters for each model (Local) and non-redudant
% list of parameters for all models (Global) and a linear Map to generate
% the local parameters from the global parameter vector.
kGlob = length(allParsStrict);
MapAllToEach = [];
MapConstraints = [];
StrictParsCell = cell(length(allParsStrict),1);
MapConstraints = cell(length(allParsConstr),1);
for i=1:nMods
    kLoc = 0;
    for j=1:length(modPar{i}.Names)
        kLoc = kLoc+1;
        switch modPar{i}.Names{j}
            case allParsStrict % strict constraint use lumped model.
                kk = find(strcmp(allParsStrict,modPar{i}.Names{j}));
                StrictParsCell(kk) = {[StrictParsCell{kk},modPar{i}.Values(j)]};
            case allParsConstr % loose constraint - free but include in regularization
                kc = find(strcmp(allParsConstr,modPar{i}.Names{j}));
                kGlob = kGlob+1; kk=kGlob;
                MapConstraints(kc) = {[MapConstraints{kc},kGlob]};
                ReducedParGuesses(kGlob,1) = modPar{i}.Values(j);
            otherwise % completely free
                kGlob = kGlob+1; kk=kGlob;
                ReducedParGuesses(kGlob,1) = modPar{i}.Values(j);
        end
        MapAllToEach{i}(kLoc,kk) = 1;       
    end
end

% Fix sizes of Maps to match total number of global parameters
for i = 1:length(MapAllToEach)
    if size(MapAllToEach{i},2)<kk
        MapAllToEach{i}(1,kk)=0;
    end
end

% Replace strict constrained parameters with mean of their guesses.
for i=1:length(allParsStrict)
    ReducedParGuesses(i,1) = mean(StrictParsCell{i});
end
ReducedParGuesses = log10(ReducedParGuesses);

% Formulate objective functions for each model/data combination.
for i=1:nMods
    OBJ{i} =  @(x)-get_fit_error(10.^[MapAllToEach{i}*x],LAS{i}.Loc_app_str,LAS{i}.vars_to_fit);
end

% Define function for constrain penalty calculation.
OBJ_Constraint = @(x)ComputeConstraintViolation(x,MapConstraints,parSimilarityConstraints);

% Combine all the objective functions. 
OverallObjective = @(x) ( getOverallObj(OBJ,x) + OBJ_Constraint(x) );

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
app = ssit.parest.updateModelSolveAndCompareToData(app);
fit_error = app.DataLoadingAndFittingTabOutputs.J_LogLk;
end

