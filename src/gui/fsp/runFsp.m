function app = runFsp(app)
% This function runs the FSP calculation after the "Run FSP" Button is
% selected within the GUI.
%% Provides an indicator when the FSP analysis is running for the user
d_prog_bar = [];%waitbar(0,'Running FSP Computation');

%% Change FSP Tolerances.
switch app.FspIntegratorTolerance.Value
    case 'Relaxed'
        relTol = 1e-2; absTol = 1e-4;
    case 'Moderate'
        relTol = 1e-3; absTol = 1e-6;
    case 'Strict'
        relTol = 1e-4; absTol = 1e-12;
end
app.SSITModel.fspOptions.fspIntegratorAbsTol = absTol;
app.SSITModel.fspOptions.fspIntegratorRelTol = relTol;

%% Set Initial Conditions and other user settings
app.SSITModel.initialCondition = eval(app.FspInitCondField.Value)';   % Pulls the inital conditions specified
Nsp = length(app.SSITModel.species);
app.SSITModel.customConstraintFuns = app.FspConstraintTable.Data(Nsp*2+1:end,1);
app.SSITModel.fspOptions.bounds = app.FspTabOutputs.bounds';
app.SSITModel.tSpan = unique(eval(app.FspPrintTimesField.Value));        % Pulls the time array from the app
app.SSITModel.fspOptions.fspTol = app.FspErrorTolField.Value;            % Pulls the specified error tolerance for the FSP Analysis
app.SSITModel.fspOptions.initApproxSS = app.initApproxSS.Value;   % Pull if using SS at initial time.
app.SSITModel.fspOptions.usePiecewiseFSP = app.FspPiecewiseCheckBox.Value;

%% Start FSPCalculations  
app.FspRunningStatus.Text = 'FSP is running...';
tic
app.SSITModel.solutionScheme = 'FSP';
[solutions,~,app.SSITModel]=app.SSITModel.solve;
app.FspRunningStatus.Text = ['FSP completed in ',num2str(toc,3),'s'];

%% Update bounds
app.FspTabOutputs.bounds = app.SSITModel.fspOptions.bounds';
app.FspTabOutputs.solutions = solutions.fsp;
delete(d_prog_bar); % Closes the Running Dialog Box
app.FspConstraintTable.Data(:,3) = num2cell(app.FspTabOutputs.bounds');   % Sets the bounds from constraints to the ones found in the FSP analysis
end
