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

%% Solve for probability sensitivities
[app.SensFspTabOutputs.solutions,constraintBounds,app.FspTabOutputs.stateSpace] = ...
    ssit.sensitivity.computeSensitivity(model,...
    parameters,...
    tspan,...
    fspTol,...
    initial_state,...
    1.0,...
    constraintFunctions,...
    constraintBounds,...
    0, 1, useMex, solutionMethod,app,...
    app.FspTabOutputs.stateSpace,...
    app.FspPiecewiseCheckBox.Value,...
    app.initApproxSS.Value,...
    app.ReactionsTabOutputs.varNames);

app.FspTabOutputs.bounds = constraintBounds';
app.FspConstraintTable.Data(:,3) = num2cell(app.FspTabOutputs.bounds');   % Sets the bounds from constraints to the ones found in the FSP analysis

close(d_prog_bar)
app.SensPlotTimeSlider.Limits = [min(tspan),max(tspan)];
app.SensPlotTimeSlider.Value = max(tspan);
updatePlotsInSensFsp(app)
end
