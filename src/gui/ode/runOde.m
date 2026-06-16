function app = runOde(app)
% This function finds the ode solution to the system of reactions

% tspan = eval(app.FspPrintTimesField.Value); % Pulls the array of times
% x0 = eval(app.FspInitCondField.Value)';   % Pulls the inital conditions specified
% inputs = app.ReactionsTabOutputs.inputs; % Pulls the input equations specified
% pars = app.ReactionsTabOutputs.parameters; % Pulls the parameters specified
% stoich_mat = app.ReactionsTabOutputs.stoichMatrix; % Pulls the stoichiometry from the reactions
% propens = app.ReactionsTabOutputs.propensities; % Pulls the propensity functions

app.SSITModel.solutionScheme = 'ode';
try
    [~,~,app.SSITModel] = app.SSITModel.solve;
    % app.FspTabOutputs.tOde = app.SSITModel.tSpan;
    % app.FspTabOutputs.odeSolutions = solnOde.ode;
catch
    disp('ODE Solution not found -- possible error with unsupported logical functions')
    app.SsaShowOdeCheckBox.Value = false;
end

% solve for the ode solutions
% [t_ode,ode_solutions] = ssit.moments.solveOde(x0, tspan, stoich_mat, propens, pars, inputs, app.SSITModel.species);
% define variables in the app for the ode solution
