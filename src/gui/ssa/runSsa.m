function [] = runSsa(app)
% This function runs the SSA analysis in the GUI after the user has
% selected the runSsa button. Through this analysis, the results are then
% pushed to the SSA tab for view by trajectory or histogram
% Gives indicator when the SSA button has been selected

if verLessThan('matlab','9.4') % This if statement allows for the script to be run in other versions than the latest version
    app.SsaRunButton.BackgroundColor = 'yellow';
    pause(2)
else
    f = app.UIFigure;
    d_prog_bar = uiprogressdlg(f,'Title','Running SSA Computation');
end

app.SSITModel.initialCondition = eval(app.SsaInitCondField.Value)';  % Pulls the data from the initial conditions field
app.SSITModel.tSpan = eval(app.PrintTimesEditField.Value);  % Pulls the time array from the app
app.SSITModel.ssaOptions.Nexp = 1;
app.SSITModel.ssaOptions.nSimsPerExpt = app.SsaNumSimField.Value/length(app.SSITModel.tSpan);  % Pulls the number of Simulations from the app
app.SSITModel.ssaOptions.useTimeVar = ~isempty(app.SSITModel.inputExpressions);
app.SSITModel.ssaOptions.signalUpdateRate = str2double(app.SsaSignalUpdateRateField.Value);
app.SSITModel.ssaOptions.useParallel = app.SsaParallelCheckBox.Value;
app.SSITModel.ssaOptions.useGPU = app.SsaGPUCheckBox.Value;
app.SSITModel.solutionScheme = 'SSA';

tic
app.SsaRunStatus.Text = 'Running...';    
[~,~,app.SSITModel] = app.SSITModel.solve;

app.SsaRunStatus.Text = ['Complete in ',num2str(toc),' seconds'];
% Re-sets the running indicator from the beginning of the function
if verLessThan('matlab','9.4') % This if statement allows for user to know it is ready to go for any version of Matlab
    app.SsaRunButton.BackgroundColor = 'green';
    pause(2)
else
    close(d_prog_bar); % Closes the Running Dialog Box
end
% app.StochasticSimulationTabOutputs.samples = app.SSITModel.Solutions.trajs; % stores the calculated Trajectory in property
% app.StochasticSimulationTabOutputs.PDOsamples = app.SSITModel.Solutions.trajsDistorted; % stores the calculated Trajectory in property
end
