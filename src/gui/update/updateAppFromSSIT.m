function updateAppFromSSIT(app,saveCopy)
arguments
    app
    saveCopy=false
end
% Update the app properties based on the SSIT data

%% GUI Information
app.ReactionsTabOutputs.stoichMatrix = app.SSITModel.stoichiometry;
if isempty(app.SSITModel.parameters)
    app.SSITModel.parameters = cell(0,2);
end
app.ReactionsTabOutputs.parameters = app.SSITModel.parameters;
app.ReactionsTabOutputs.inputs = app.SSITModel.inputExpressions;
app.ReactionsTabOutputs.propensities = app.SSITModel.propensityFunctions;
app.ReactionsTabOutputs.presetParameters = [app.SSITModel.parameters{:,2}];
if isempty(app.SSITModel.inputExpressions)
    app.SSITModel.inputExpressions = cell(0,2);
end
app.ReactionsTabOutputs.presetInputs = app.SSITModel.inputExpressions(:,2);

%% Reactions
app.ModelReactionTable.Data = {};
nRxn = size(app.SSITModel.stoichiometry,2);
nSp =size(app.SSITModel.stoichiometry,1);
for iRxn = 1:nRxn
    reactants = '';
    products = '';
    for jSp = 1:nSp
        if app.SSITModel.stoichiometry(jSp, iRxn) < 0
            reactants = append(reactants, app.SSITModel.species{jSp},'(',num2str(-app.SSITModel.stoichiometry(jSp, iRxn)),'),');
        elseif app.SSITModel.stoichiometry(jSp, iRxn) > 0
            products = append(products, app.SSITModel.species{jSp},'(',num2str(app.SSITModel.stoichiometry(jSp, iRxn)),'),');
        end
    end
    if isempty(reactants); reactants = '-,';end
    if isempty(products); products = '-,';end   
    app.ModelReactionTable.Data{iRxn, 1} = append('R',num2str(iRxn));
    app.ModelReactionTable.Data{iRxn, 2} = reactants(1:end-1);
    app.ModelReactionTable.Data{iRxn, 3} = products(1:end-1);
    app.ModelReactionTable.Data(iRxn, 4) = app.SSITModel.propensityFunctions(iRxn);
end

app.DeleteReactionDropDown.Items = {'Select'};
for iRxn = 1:nRxn
    app.DeleteReactionDropDown.Items{iRxn+1} = ['R',num2str(iRxn)];
end
app.DeleteReactionDropDown.Value=app.DeleteReactionDropDown.Items{1};
%% Parameter Table
app.ModelParameterTable.Data = app.SSITModel.parameters;
app.SensParameterSelectionTable.Data = app.SSITModel.parameters;
app.fit_parameters_table.Data ={};
for iPar = 1:size(app.SensParameterSelectionTable.Data,1)
    app.SensParameterSelectionTable.Data{iPar,3} = 'y';
    app.fit_parameters_table.Data(iPar,1) = app.SSITModel.parameters(iPar,1);
    app.fit_parameters_table.Data(iPar,2) = app.SSITModel.parameters(iPar,2);
end
app.fit_parameters_table.Data(:,3) = {'n'};
app.fit_parameters_table.Data(app.SSITModel.fittingOptions.modelVarsToFit,3) = {'y'};

% FIM Options
FIMMetrics = {'Determinant','Smallest Eigenvalue','Trace'};
FIMMetrics(end+1:end+size(app.SSITModel.parameters,1)) = app.SSITModel.parameters(:,1);
app.FIMMetricorParameterDropDown.Items = FIMMetrics;
app.FIMMetricorParameterDropDown.Value = 'Determinant';

%% Input Expressions
nSig = size(app.SSITModel.inputExpressions,1);
app.ModelInputTable.Data = {};
for iSig = 1:nSig
    app.ModelInputTable.Data{iSig,1} = app.SSITModel.inputExpressions{iSig,1};
    app.ModelInputTable.Data{iSig,2} = app.SSITModel.inputExpressions{iSig,2};
end

app.DeleteInputDropDown.Items = {'Select'};
if nSig>=1
    app.DeleteInputDropDown.Enable=true;
    for iSig = 1:nSig
        app.DeleteInputDropDown.Items{iSig+1} = app.SSITModel.inputExpressions{iSig,1};
    end
else
    app.DeleteInputDropDown.Enable=false;
end
app.DeleteInputDropDown.Value=app.DeleteInputDropDown.Items{1};

%% Information
app.ModelAbout.Value = app.SSITModel.description;

%% Initial Conditions
app.SsaInitCondField.Value = append('[',num2str(app.SSITModel.initialCondition'),']');
app.FspInitCondField.Value = append('[',num2str(app.SSITModel.initialCondition'),']');
app.SensInitCondField.Value = append('[',num2str(app.SSITModel.initialCondition'),']');

%% Model Solution Settings
app.PrintTimesEditField.Value = append('[',num2str(app.SSITModel.tSpan),']');
app.FspPrintTimesField.Value = append('[',num2str(app.SSITModel.tSpan),']');
app.SensPrintTimesEditField.Value = append('[',num2str(app.SSITModel.tSpan),']');
app.ListofMeasurementTimesEditField.Value = append('[',num2str(app.SSITModel.tSpan),']');

%% Update user dropdown choices
speciesLists = {'SpeciestoShowListBox_2','SpeciestoShowListBoxMargFSP',...
    'SpeciestoShowListBoxMargFSPvT','SpeciestoShowListBoxMeans','JointSp1','JointSp2',...
    'SpeciestoShowListBox','SpeciesForSensPlot','ObservableSpeciesListBox',...
    'SpeciesForFitPlot'};
for sp = speciesLists
    app.(sp{1}).Items = app.SSITModel.species;
end
app.SpeciesTable.Data ={};
for iSp = 1:size(app.SSITModel.species,1)
    app.SpeciesTable.Data{iSp,1} = app.SSITModel.species{iSp};
    app.SpeciesTable.Data{iSp,2} = app.SSITModel.initialCondition(iSp);
end   

app.DeleteSpeciesDropDown.Items = {'Select'};
for iSp = 1:size(app.SSITModel.species,1)
    app.DeleteSpeciesDropDown.Items{iSp+1} = app.SSITModel.species{iSp};
end
app.DeleteSpeciesDropDown.Value=app.DeleteSpeciesDropDown.Items{1};


%% FSP Options
makeDefaultConstraints(app);
readConstraintsForAdaptiveFsp(app);

%% Add data elements
if ~isempty(app.SSITModel.dataSet)
    for i = 1:length(app.SSITModel.dataSet.times)
        app.ParEstFitTimesList.Items{i} = num2str(app.SSITModel.dataSet.times(i));
    end
end

%% Update app model description
summarizeModelReactions(app.SpeciesTable.Data(:,1),app.ModelReactionTable.Data,app.ModelInputTable.Data,app.ReactionTextAxes)

%% Update propensity function codes
app.ModelhasnotbeenupdatedLabel.Text = 'Writing Propensity Function Codes.';  % Fills in update for user
drawnow

app.SSITModel = app.SSITModel.formPropensitiesGeneral(app.ModelFile.modelName);
app.ModelhasnotbeenupdatedLabel.Text = 'Model is up to date.';  % Fills in update for user

if saveCopy
    fileName = app.ModelFile.fileName;
    if strcmp(fileName(end-1:end),'.m')
        fileName = append(fileName(1:end-2),'.mat');
    elseif strcmp(fileName(end-3:end),'.mat')
        fileName = append(fileName(1:end-4),'_SSIT.mat');
    end
    eval(append(app.ModelFile.modelName,' = app.SSITModel;'));
    save(fileName,app.ModelFile.modelName)
end



end
