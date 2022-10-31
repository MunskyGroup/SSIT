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

model = app.Model;
S = model.stoichiometry;
inputs = model.timeVaryingInputExpressions;
pars = app.ReactionsTabOutputs.parameters;
W = model.processPropensityStrings(app.ReactionsTabOutputs.propensities,inputs,pars,'fun',app.ReactionsTabOutputs.varNames);

x0 = eval(app.SsaInitCondField.Value)';  % Pulls the data from the initial conditions field
timesToRecord = eval(app.PrintTimesEditField.Value);  % Pulls the time array from the app

Nt = length(timesToRecord);
Nsim = app.SsaNumSimField.Value;  % Pulls the number of Simulations from the app
nSpecies = length(app.ReactionsTabOutputs.varNames);
Trajs = zeros(nSpecies,Nt,Nsim);                       % Creates an empty Trajectories matrix from the size of the time array and number of simulations
useTimeVar = ~isempty(app.ReactionsTabOutputs.inputs);
signalUpdateRate = str2double(app.SsaSignalUpdateRateField.Value);

tic
app.SsaRunStatus.Text = 'Running...';    

if app.SsaParallelCheckBox.Value
    try
        parpool
    catch
    end
    parfor isim = 1:Nsim
        Trajs(:,:,isim) = ssit.ssa.runSingleSsa(x0,S,W,timesToRecord,useTimeVar,signalUpdateRate);
    end
else
    for isim = 1:Nsim
        Trajs(:,:,isim) = ssit.ssa.runSingleSsa(x0,S,W,timesToRecord,useTimeVar,signalUpdateRate);
    end
end
app.SsaRunStatus.Text = ['Complete in ',num2str(toc),' seconds'];
% Re-sets the running indicator from the beginning of the function
if verLessThan('matlab','9.4') % This if statement allows for user to know it is ready to go for any version of Matlab
    app.SsaRunButton.BackgroundColor = 'green';
    pause(2)
else
    close(d_prog_bar); % Closes the Running Dialog Box
end
app.StochasticSimulationTabOutputs.samples = Trajs; % stores the calculated Trajectory in property
end
