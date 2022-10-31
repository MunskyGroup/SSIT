function [t_ode,ode_solutions] = solveOde(x0, tspan, stoichMatrix, propens, pars, inputs, varNames)
% ode_solutions = solve_ode(x0, tspan, stoichMatrix, propens, pars,
% inputs) solves the deterministic ODE model.
% x0: initial state. This must be a column vector.
% tspan: array of output times.
% stoichMatrix: stoichiometry matrix, each column represents the net change
% of one reaction.
% propstrings: cells of propensity strings, using the same format as in the
% Reaction tab of the app
% pars: 2D cell of parameters, first column contain the names of
% parameters, second column contains the values.
% inputs: cell of time-varying inputs

propstrings = ssit.SrnModel.processPropensityStrings(propens,inputs,pars,"str",varNames);
n_reactions = length(propstrings);
propensities = cell(n_reactions, 1);
for i = 1:n_reactions
    propstrings{i} = strrep(propstrings{i},'..','.');
    propensities{i} = ssit.Propensity.createFromString(propstrings{i}, stoichMatrix(:,i), i);
end

%% Define the right hand side of the ODE model
ode_rhs = @(t, x) stoichMatrix*generate_propensity_vector(t, x, propensities);

%% Return your output to ode_solutions
[t_ode,ode_solutions] = ode45(ode_rhs,tspan,x0);
end

function y = generate_propensity_vector(t, x, propensities)
% Helper function to output a column vector of propensity values evaluated at the
% state x at time t.

y = zeros(length(propensities), 1);
for i = 1:length(propensities)    
%     y(i) = propensities{i}.eval(t,x);
    if propensities{i}.isFactorizable
        y(i) = propensities{i}.timeDependentFactor(t).*propensities{i}.stateDependentFactor(x);
    else
        y(i) = propensities{i}.jointDependentFactor(t,x);
    end
end
end