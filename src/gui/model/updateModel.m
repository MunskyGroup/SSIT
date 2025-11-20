function updateModel(app,saveCopy,fileName,useOldXspeciesVersion)
arguments
    app
    saveCopy = false
    fileName = []
    useOldXspeciesVersion=false
end
% This function changes updates the propensity vector, parameter table, and
% inputs of the model within the GUI. This occurs after the
% UpdateModelButton had been pushed.
% TODO -- this is an outdated code that will no longer be needed once we
% move to the SSIT command line approach. However, we need to update old
% models using this to parse them into the new format.  See end of code for
% hopw this might be done.


Mod = app.ModelReactionTable.Data;    % Assigns a variable to the model set in the Reactions Table
if size(Mod,2)>4
    Mod = Mod(:,1:4);
end
% J = [];                     % Pre-allocates an empy variable to be updated in the for-loop
% for i=1:size(Mod,1)         % For-loop to build a vector of indices
%     if Mod{i,5}~='y'        % If-statement that excludes all rows that are noted as removed: 'y' indicates that it needs to be excluded from the model
%         J = [J,i];          % This builds the vector
%     end
% end
% Mod = Mod(J,:);             % Uses the index vector, J, to only include the rows that are not being removed
% app.ModelReactionTable.Data = Mod;    % Re-assigns the updated model table by removing the unwanted rows
% pars = {'x1';'x2';'x3'};    % Creates cells for parameters
inputs = {};                % Creates an empty cell for input vector
   
nRxn = size(app.ModelReactionTable.Data,1);        % Gives the number of rows in the Reactions Table
% Add species names to pars list.

if useOldXspeciesVersion
    species = {};
    for iRxn=1:nRxn                                       % For-loop over all rows in the Reactions Table
        for col = 2:3
            Rcts = Mod{iRxn,col};                               % Assigns a variable to represent the value of row i of the Reactant in the Reactions Table
            J = strfind(Rcts,'x');                           % Returns a vector of the starting indices of any occurrences of x in the Reactants variable
            J2 = strfind(Rcts,'(')-1;                        % Returns a vector of the starting indices of any occurrences of x in the Reactants variable
            for j=1:length(J)                                % For-loop over the number of reactants found by J
                species(end+1,1) = {Rcts(J(j):J2(j))};
            end
        end
    end
    species = unique(species,'stable');
    species(:,2) = num2cell(zeros(length(species),1));
else
    species = app.SSITModel.species;
    % app.SSITModel.initialCondition = [app.SpeciesTable.Data{:,2}]
    % species(:,2) = app.SpeciesTable.Data{:,2};%num2cell(app.SSITModel.initialCondition);    
end
% pars = unique(species,'stable');    % Removes any repeated propensity
% app.ReactionsTabOutputs.varNames = pars;
% app.SpeciestoShowListBox.Items = pars;
% app.SpeciestoShowListBox_2.Items = pars;
% app.SpeciestoShowListBoxMargFSP.Items = pars;
% app.SpeciestoShowListBoxMargFSPvT.Items = pars;
% app.SpeciestoShowListBoxMeans.Items = pars;
% app.ObservableSpeciesListBox.Items = pars;
% app.ObservableSpeciesListBox.Value = pars;
% app.SpeciesForFitPlot.Items = pars;
% app.SpeciesForSensPlot.Items = pars;
% app.JointSp1.Items = pars;
% if length(pars)>1
%     app.JointSp2.Enable = 1;
%     app.PlotJointDistributionsButton.Enable = 1;
%     app.JointSp2.Items = pars; app.JointSp2.Value=pars{2};
% else
%     app.JointSp2.Enable = 0;
%     app.PlotJointDistributionsButton.Enable = 0;
% end

nSpecies = size(species,1);
pars = species(:,1);
if size(Mod,1)>=1                                   % If-statement to ensure that the table is not empty    
    S = zeros(nSpecies,nRxn);                      % Assigns a zeros ReactionsTabOutputs.stoichMatrixiometry matrix for the three species to be filled in with the following for-loops
    for iRxn=1:nRxn                                       % For-loop over all rows in the Reactions Table
        %% Build ReactionsTabOutputs.stoichMatrixiometry matrix
        Rcts = Mod{iRxn,2};                                 
        Prods = Mod{iRxn,3};
        for iSp = 1:nSpecies
            % Reactants
            J = strfind(Rcts,species{iSp,1});
            J1 = strfind(Rcts,'(');
            J2 = strfind(Rcts,')');
            for j=1:length(J)
                S(iSp,iRxn)=S(iSp,iRxn)-str2num(Rcts(J1(j)+1:J2(j)-1));      % Assign reaction effect on stoichiometry
            end

            % Products
            J = strfind(Prods,species{iSp,1});
            J1 = strfind(Prods,'(');
            J2 = strfind(Prods,')');
            for j=1:length(J)                                % For-loop over the number of products found by J
                S(iSp,iRxn)=S(iSp,iRxn)+str2num(Prods(J1(j)+1:J2(j)-1)); % Assign product effect on stoichiometry
            end
        end
        
        %% Define Propensity Functions
        % Propensity Function Parameters
        pars = [pars(:);symvar(Mod{iRxn,4})];
        pars = unique(pars,'stable');    % Removes any repeated parameter
        
        % Change Propensity Functions to allow for Vector ReactionsTabOutputs.inputs
        op_var = {'*','/','^'};
        Props = Mod{iRxn,4};                                   % Assigns a variable to represent the value of reaction i of the propensity function - %% May be reduntant and can remove. Already degined on line 214 and 223
        for jop = 1:length(op_var)
            op = op_var{jop};
            Props = strrep(Props, op, ['.' op]);
            while contains(Props,'..')
                Props = strrep(Props, '..', '.');
            end
        end
        Props_vec{iRxn}=Props;                                 % Re-defines the propensity vector based on the results of the two for-loops above
    end
    
    %% Check Input Signals for additional parameters
    for iRxn = 1:size(app.ModelInputTable.Data,1)
        pars = [pars(:);symvar(app.ModelInputTable.Data{iRxn,2})];
        pars = unique(pars,'stable');    % Removes any repeated propensity
    end
    J = ~strcmp(pars,'t');
    pars = pars(J);
    
    %% Sort Parameters Into Constants and ReactionsTabOutputs.inputs
    pars = pars(nSpecies+1:end);
    
    % Re-order parameters from longest to shortest:
    [~,reorderPars] = sort(cellfun(@length, pars),'descend');
    pars = pars(reorderPars);
    
    J = zeros(length(pars),1,'logical');
    if size(app.ModelInputTable.Data,1)
        [~,Jn] = intersect(pars,app.ModelInputTable.Data{:,1});
        J(Jn) = 1;
        if sum(J)<size(app.ModelInputTable.Data,1)
            disp('WARNING - Removing one or more unused signals.')
        end
    end
    % J = startsWith(pars,'I');
    inputs = pars(J);
    pars = pars(~J);
    
    %% Define propensity functions as inline statements.
    % Create list of parameters
    tmp =[];
    if length(pars)>=1
        tmp = pars{1};
        for j=2:length(pars)
            tmp = [tmp,',',pars{j}];
        end
    end
    % Create list of inputs
    tmp2 =[];
    if length(inputs)>=1
        tmp2 = [inputs{1}];
        for j=2:length(inputs)
            tmp2 = [tmp2,',',inputs{j}];
        end
    end
    
    % Concatenate lists
    if isempty(tmp)
        tmp = tmp2;
    elseif ~isempty(tmp2)
        tmp = [tmp,',',tmp2];
    end
    
    % This section loads in the input signals and adds relevant parameters
    % to the list.
    if isempty(inputs)
        inputs=zeros(0,2);
        app.ReactionsTabOutputs.inputs=inputs;
        app.ReactionsTabOutputs.presetInputs=[];
    else
        if ~isempty(app.ReactionsTabOutputs.inputs)                                % Logical command to execute if there was previously defined a set of inputs
            if isempty(app.ReactionsTabOutputs.presetInputs)                                % Logical command if there were no pre-defined parameters from the model
                inputs(:,2) = {'1+sin(t)'};                         % Fills the input vector with the default time dependence character vector
            else
                [~,iold,inew] = intersect(app.ReactionsTabOutputs.inputs(:,1),inputs);  % Returns index vectors of where the values of the new inputs overlap with the old ones
                inputs(:,2) = {'1+sin(t)'};                         % Fills the input vector with the default time dependence character vector
                inputs(inew,2) = app.ReactionsTabOutputs.presetInputs(iold);               % Replaces the ones in the new values with the previously defined parameters
            end

        else                                                  % Defines the input vector when there was no previous run
            if ~isempty(app.ReactionsTabOutputs.presetInputs)                               % Logical command that checks if there are pre-defined inputs from the preset models
                inputs(:,2) = app.ReactionsTabOutputs.presetInputs;                         % Re-defines the input values if there was a pre-defined set of inputs from the preset models
            else
                inputs(:,2) = {'1+sin(t)'};                             % Defines the inputs vector as the default if there was no pre-defined set of parameters from the preset models
            end
        end
    end
    
    for iRxn = 1:size(inputs,1)
        pars = [pars(:);symvar(inputs{iRxn,2})];
        pars = unique(pars,'stable');    % Removes any repeated propensity
    end
    J = ~strcmp(pars,'t');
    pars = pars(J);
        
    % The following corrects for the removing of parameters in the interface
    if ~isempty(app.ReactionsTabOutputs.parameters)                              % This checks if the parameters vector is not empty
        if isempty(app.ReactionsTabOutputs.presetParameters)                       % Logical command to execute if the previous parameter values are empty
            pars(:,2) = {1};                           % Fills the empty parameters vector with ones
        else
            [~,iold,inew] = intersect(app.ReactionsTabOutputs.parameters(:,1),pars);  % Returns index vectors of where the values of the new parameters overlap with the old ones
            pars(:,2) = {1};                               % Fills all parameters vector with ones
            if iscell(app.ReactionsTabOutputs.presetParameters(iold))
                pars(inew,2) = app.ReactionsTabOutputs.presetParameters(iold);               % Replaces the ones in the new values with the previously defined parameters
            else
                pars(inew,2) = num2cell(app.ReactionsTabOutputs.presetParameters(iold));               % Replaces the ones in the new values with the previously defined parameters
            end
        end
    else                                               % This defines the parameter vector when there was none before
        if ~isempty(app.ReactionsTabOutputs.presetParameters)                       % This if statement checks if there are pre-defined parameters from the preset examples
            pars(:,2) = app.ReactionsTabOutputs.presetParameters;                       % Re-defines the parameter values if there was a pre-defined set of parameters from the preset models
        else
            pars(:,2) = {1};                                % Defines the parameters vector as ones if there was no pre-defined set of parameters from the preset models
        end
    end
    
    % if isempty(app.ReactionsTabOutputs.citations)                % Identifies if there is citation for a preset example
    %     app.ModelCitationforPresetExamplesTextArea.Value = 'Citation for Preset Example Unavailable';                       % Sets the Citation as Unavailable if none was specified
    % else
    %     app.ModelCitationforPresetExamplesTextArea.Value = sprintf('Citation for Preset Example: %s\n',char(app.ReactionsTabOutputs.citations)); % Sets the Citation as the given citation from the preset example
    % end
    
%     if ~isempty(app.ReactionsTabOutputs.initialCondition)           % Identifies the initial conditions given from the preset example
%         app.SsaInitCondField.Value = app.ReactionsTabOutputs.initialCondition;       % Sets the initial conditions in the SSA Tab
%         app.FspInitCondField.Value = app.ReactionsTabOutputs.initialCondition;     % Sets the initial conditions in the FSP Tab
%     else
        app.SsaInitCondField.Value = ['[',num2str(zeros(1,nSpecies)),']'];                % Sets the initial conditions in the SSA Tab to default if unspecified in the model
        app.FspInitCondField.Value = ['[',num2str(zeros(1,nSpecies)),']'];              % Sets the initial conditions in the SSA Tab to default if unspecified in the model
        app.SensInitCondField.Value = ['[',num2str(zeros(1,nSpecies)),']'];              % Sets the initial conditions in the SSA Tab to default if unspecified in the model
%     end
    
    % app.ModelInputTable.Data = inputs;       % Fills in the Input Table
    % app.ModelParameterTable.Data = pars;                % Fills in the Parameters Table
end

app.ReactionsTabOutputs.stoichMatrix = S;
app.ReactionsTabOutputs.parameters = pars;
app.ReactionsTabOutputs.inputs = inputs;
app.ReactionsTabOutputs.propensities = Props_vec;
app.ReactionsTabOutputs.presetParameters = [pars{:,2}];
app.ReactionsTabOutputs.presetInputs = inputs(:,2);
% selectSensitivityParametersTable(app);                          % Updates the Sensitivity Analysis table with the Parameters from the Reactions tab
% editSelectParametersTable(app);                                 % Creates the Selectivity Analysis Vectors

% app.Model = ssit.SrnModel(S, Props_vec, pars(:,1), inputs);

% Clears the graphs in the SSA and FSP Tabs
cla(app.SsaTrajAxes);
cla(app.SsaHistAxes);
cla(app.FspAxes);
% cla(app.MomentsAxes);
% cla(app.MomentsHistAxes);

if ~isempty(app.ReactionsTabOutputs.presetInputs)
    app.SsaSignalUpdateRateField.Enable = 1;
else
    app.SsaSignalUpdateRateField.Enable = 0;
end

%% Reset the Propensity Functions
% if isfield(app.FspTabOutputs,'PropensitiesGeneral')
%     app.FspTabOutputs = rmfield(app.FspTabOutputs,'PropensitiesGeneral');
% end

%% Reset the sensitivity and FIM calculations for new model
app.SensFspTabOutputs.solutions = [];
% FIMMetrics = {'Determinant','Smallest Eigenvalue','Trace'};
% FIMMetrics(end+1:end+size(pars,1)) = pars(:,1);
% app.FIMMetricorParameterDropDown.Items = FIMMetrics;
% app.FIMMetricorParameterDropDown.Value = 'Determinant';
% app.FIMSpecies1.Items = pars(:,1)'; app.FIMSpecies1.Value = pars{1,1};
% app.FIMSpecies2.Items = pars(:,1)'; 
% if size(pars,1)>=2
%     app.FIMSpecies2.Value = pars(2,1);
% else
%     app.FIMSpecies2.Value = pars(1,1);
% end
% app.ModelUncertaintyDropDown.Value = 'Use Fixed Model';

app.FIMTabOutputs.FIMPrior = [];
app.FIMTabOutputs.CellsPerTimePoint = [];
app.FspTabOutputs.odeSolutions = [];
% makeDefaultConstraints(app);
% readConstraintsForAdaptiveFsp(app);

%% Update table of species names:
% app.NameTable.Data={};
% app.NameTable.Data(:,1) = app.ReactionsTabOutputs.varNames;
% app.NameTable.Data(:,2) = app.ReactionsTabOutputs.varNames;

%% Update SSIT model
app.SSITModel = SSIT('Empty');
app.SSITModel.species = species(:,1);
app.SSITModel.stoichiometry = S;
app.SSITModel.parameters = pars;
app.SSITModel.propensityFunctions = Props_vec';
app.SSITModel.inputExpressions = inputs;
% if ~isempty(species)
    app.SSITModel.initialCondition = [species{:,2}]';
% end
app.SSITModel.description = app.ModelAbout.Value;
% k = strfind(fileName,'.'); k=k(end);
% k1 = strfind(fileName,'/'); k1=k1(end);
% ModelName = fileName(k1+1:k-1);
% disp('Writing Propensity Function Codes.')
% app.SSITModel = app.SSITModel.formPropensitiesGeneral(app.ModelFile.modelName);
updateAppFromSSIT(app,true);

