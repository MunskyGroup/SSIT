function app = runSensitivity(app,useProgBar)
%RUNSENSITIVITY
arguments
    app
    useProgBar = true;
end
if useProgBar
    try
        f = app.UIFigure;
        d_prog_bar = uiprogressdlg(f,'Title','Running Sensitivity Computation','Indeterminate','on');
    catch
        d_prog_bar=[];
    end
else
    d_prog_bar=[];
end
app.SensParDropDown.Items = app.ReactionsTabOutputs.parameters(:,1);
app.SensParDropDown.ItemsData = app.ReactionsTabOutputs.parameters(:,1);
app.SensParDropDown.Value = app.ReactionsTabOutputs.parameters{1,1};

if app.SensitivityFunctionButton.Value == 1
    solutionMethod = 'forward';
else
    solutionMethod = 'finitediff';
end

useMex = 0;%app.SensUseSundialCheckBox.Value;
%% Read parameter values from app
model = app.Model;
initial_state = eval(app.FspInitCondField.Value)';
tspan = unique(eval(app.SensPrintTimesEditField.Value));        % Pulls the time array from the app
fspTol = app.SensErrorTolEditField.Value;            % Pulls the specified error tolerance for the FSP Analysis
if isempty(app.FspTabOutputs.bounds) % if we have yet set up the bounds in FSP tab, resort to default bounds
    makeDefaultConstraints(app);
    readConstraintsForAdaptiveFsp(app);
end
constraintFunctions = readConstraintsForAdaptiveFsp(app);
constraintBounds = app.FspTabOutputs.bounds';

parameters = app.SensParameterSelectionTable.Data(:,1:2);

n_reactions = size(app.ReactionsTabOutputs.stoichMatrix,2);

if isfield(app.FspTabOutputs,'PropensitiesGeneral')&&~isempty(app.FspTabOutputs.PropensitiesGeneral{1}.sensTimeFactorVec)
    PropensitiesGeneral = app.FspTabOutputs.PropensitiesGeneral;
else
    sm = cell(1,n_reactions);
    logicTerms = cell(1,n_reactions);
    logCounter = 0;
    for i = 1:n_reactions
        st = app.ReactionsTabOutputs.propensities{i};
        for jI = 1:size(app.ReactionsTabOutputs.inputs,1)
            st = strrep(st,app.ReactionsTabOutputs.inputs{jI,1},['(',app.ReactionsTabOutputs.inputs{jI,2},')']);
        end
        [st,logicTerms{i},logCounter] = ssit.Propensity.stripLogicals(st,app.ReactionsTabOutputs.varNames,logCounter);
        sm{i} = str2sym(st);
    end

    % if obj.useHybrid
    %     PropensitiesGeneral = ssit.Propensity.createAsHybridVec(sm, obj.stoichiometry,...
    %         obj.parameters, obj.species, obj.hybridOptions.upstreamODEs, logicTerms, prefixName);
    % else
    PropensitiesGeneral = ssit.Propensity.createAsHybridVec(sm, ...
        app.ReactionsTabOutputs.stoichMatrix,...
        app.ReactionsTabOutputs.parameters, ...
        app.ReactionsTabOutputs.varNames, ...
        [], logicTerms, 'TEMPPropensSens',true);
    app.FspTabOutputs.PropensitiesGeneral = PropensitiesGeneral;
end

if isfield(app.FspTabOutputs,'solutions')
    fspSoln.fsp = app.FspTabOutputs.solutions;
    fspSoln.stateSpace = app.FspTabOutputs.stateSpace;
else
    fspSoln.fsp = [];
end

%% Solve for probability sensitivities
[app.SensFspTabOutputs.solutions,constraintBounds,app.FspTabOutputs.stateSpace] = ...
    ssit.sensitivity.computeSensitivity(...
    parameters,...
    PropensitiesGeneral,...
    tspan,...
    fspTol,...
    initial_state,...
    1.0,...
    app.ReactionsTabOutputs.stoichMatrix,...
    constraintFunctions,...
    constraintBounds,...
    0, 1, useMex, solutionMethod,app,...
    app.FspTabOutputs.stateSpace,...
    app.FspPiecewiseCheckBox.Value,...
    app.initApproxSS.Value,...
    app.ReactionsTabOutputs.varNames,...
    true,fspSoln,false,[],false,[],[],[]);

app.FspTabOutputs.bounds = constraintBounds';
app.FspConstraintTable.Data(:,3) = num2cell(app.FspTabOutputs.bounds');   % Sets the bounds from constraints to the ones found in the FSP analysis

close(d_prog_bar)
app.SensPlotTimeSlider.Limits = [min(tspan),max(tspan)];
app.SensPlotTimeSlider.Value = max(tspan);
updatePlotsInSensFsp(app)
end
