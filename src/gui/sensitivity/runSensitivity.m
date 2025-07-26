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
    app.SSITModel.sensOptions.solutionMethod = 'forward';
else
    app.SSITModel.sensOptions.solutionMethod = 'finitediff';
end

%% Set Initial Conditions and other user settings
app.SSITModel.initialCondition = eval(app.FspInitCondField.Value)';   % Pulls the inital conditions specified
Nsp = length(app.SSITModel.species);
app.SSITModel.customConstraintFuns = app.FspConstraintTable.Data(Nsp*2+1:end,1);
app.SSITModel.fspOptions.bounds = app.FspTabOutputs.bounds';
app.SSITModel.tSpan = unique(eval(app.FspPrintTimesField.Value));        % Pulls the time array from the app
app.SSITModel.fspOptions.fspTol = app.FspErrorTolField.Value;            % Pulls the specified error tolerance for the FSP Analysis

%% Solve for probability sensitivities
app.SSITModel.solutionScheme = 'fspsens';
app.FspRunningStatus.Text = 'FSP-Sens is running...';
tic
[solutions,constraintBounds] = app.SSITModel.solve;
app.FspRunningStatus.Text = ['FSP-Sens completed in ',num2str(toc,3),'s'];

app.FspTabOutputs.bounds = constraintBounds';
app.SensFspTabOutputs.solutions = solutions.sens;
app.FspConstraintTable.Data(:,3) = num2cell(app.FspTabOutputs.bounds');   % Sets the bounds from constraints to the ones found in the FSP analysis

close(d_prog_bar)
app.SensPlotTimeSlider.Limits = [min(app.SSITModel.tSpan),max(app.SSITModel.tSpan)];
app.SensPlotTimeSlider.Value = max(app.SSITModel.tSpan);
updatePlotsInSensFsp(app)
end
