function updateAppFromSSIT(app)
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

%% Parameter Table
app.ModelParameterTable.Data = app.SSITModel.parameters;

%% Input Expressions
nSig = size(app.SSITModel.inputExpressions,1);
app.ModelInputTable.Data = {};
for iSig = 1:nSig
    app.ModelInputTable.Data{iSig,1} = app.SSITModel.inputExpressions{iSig,1};
    app.ModelInputTable.Data{iSig,2} = app.SSITModel.inputExpressions{iSig,2};
end

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

end
