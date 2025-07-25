function app = runOde(app)
% This function finds the ode solution to the system of reactions
% TODO - This function currently uses an older version of the ODE solver.
% It works, but should be replaced with the newer one that is built into
% the SSIT.  

tspan = eval(app.FspPrintTimesField.Value); % Pulls the array of times
x0 = eval(app.FspInitCondField.Value)';   % Pulls the inital conditions specified
inputs = app.ReactionsTabOutputs.inputs; % Pulls the input equations specified
pars = app.ReactionsTabOutputs.parameters; % Pulls the parameters specified
stoich_mat = app.ReactionsTabOutputs.stoichMatrix; % Pulls the stoichiometry from the reactions
propens = app.ReactionsTabOutputs.propensities; % Pulls the propensity functions

% solve for the ode solutions
[t_ode,ode_solutions] = ssit.moments.solveOde(x0, tspan, stoich_mat, propens, pars, inputs, app.SSITModel.species);
% define variables in the app for the ode solution
app.FspTabOutputs.tOde = t_ode;
app.FspTabOutputs.odeSolutions = ode_solutions;