function [t_ode,ode_solutions] = solveOde2(x0, tspan, stoichMatrix, propensities, useSSIC)
% ode_solutions = solve_ode(x0, tspan, stoichMatrix, propens, pars,
% inputs) solves the deterministic ODE model.
% x0: initial state. This must be a column vector.
% tspan: array of output times.
% stoichMatrix: stoichiometry matrix, each column represents the net change
% of one reaction.

% propstrings = ssit.SrnModel.processPropensityStrings(propens,inputs,pars,"str",varNames);
% n_reactions = length(propens);
% propensities = cell(n_reactions, 1);
% for i = 1:n_reactions
%     propstrings{i} = strrep(propstrings{i},'..','.');
%     propensities{i} = ssit.Propensity.createFromString(propstrings{i}, stoichMatrix(:,i), i);
% end
arguments
    x0
    tspan
    stoichMatrix
    propensities
    useSSIC = false
end

%% Define the right hand side of the ODE model
ode_rhs = @(t, x) stoichMatrix*generate_propensity_vector(t, x, propensities);

if useSSIC
     OPTIONS = optimoptions('fsolve','display','none','OptimalityTolerance',1e-8,'MaxIterations',2000);
     FUN = @(x)ode_rhs(0,x);
     x0b = fsolve(FUN,x0,OPTIONS);
%     if min(x0b)<0
        FUN = @(t,x)ode_rhs(0,x);
        [~,ode_solutions] = ode45(FUN,max(tspan)*[0,500,1000],x0b);
        x0b = ode_solutions(end,:);
%     else
%         b = pinv(stoichMatrix)*(x0b-x0);
%         if max(abs(stoichMatrix*b-(x0b-x0)))>=1e-10
%             FUN = @(~,x)ode_rhs(0,x);
%             [~,ode_solutions] = ode45(FUN,max(tspan)*[0,500,1000],x0);
%             x0b = ode_solutions(end,:);
%         end
%     end
    x0 = x0b;
end

%% Return your output to ode_solutions
[t_ode,ode_solutions] = ode15s(ode_rhs,tspan,x0);
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