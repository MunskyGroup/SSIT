function [t_ode,ode_solutions] = solveOde2(x0, tspan, stoichMatrix, propensities, parameters, useSSIC)
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
    parameters
    useSSIC = false
end

%% Define the right hand side of the ODE model
ode_rhs = @(t,x)stoichMatrix*propensities{1}.hybridFactorVector(t,parameters,x');
% ode_rhs = @(t, x) stoichMatrix*generate_propensity_vector(t, x, propensities, parameters);
% ode_jac = @(t,x)stoichMatrix*propensities{1}.DhybridFactorDodesVec(t,parameters,x');
% ode_jac = @(t,x)stoichMatrix*generate_propensity_jacobian(t, x, propensities, parameters);

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
maxstep = min(tspan(2:end)-tspan(1:end-1))/2;
options = odeset(RelTol=1e-6,AbsTol=1e-10,MaxStep=maxstep);
[t_ode,ode_solutions] = ode45(ode_rhs,tspan,x0,options);

end

% function w = generate_propensity_vector(t, x, propensities, parameters)
% Helper function to output a column vector of propensity values evaluated at the
% state x at time t.

% w = zeros(length(propensities), 1);
% if propensities{1}.isFactorizable
    % w = propensities{1}.hybridFactorVector(t,parameters,x');
    % wx = propensities{1}.xFactorVector(x,parameters);
    % w = wt.*wx;
% else
    % for i = 1:length(propensities)
        % if propensities{i}.isFactorizable
        % w(i) = wt(i).*propensities{i}.stateDependentFactor(x,parameters);
        % else
        % w(i) = propensities{i}.jointDependentFactor(t,x);
        % end
    % end
% end
% end

function dwdx = generate_propensity_jacobian(t, x, propensities, parameters)
% Helper function to output a column vector of propensity values evaluated at the
% state x at time t.

dwdx = zeros(length(propensities), size(x,1));
if propensities{1}.isFactorizable
    wt = propensities{1}.hybridFactorVector(t,parameters);
end
for i = 1:length(propensities)    
    if propensities{i}.isFactorizable
        dwdx(i,:) = wt(i).*propensities{i}.DstateFactorDstate(x,parameters);
    else
        error('Missing Code')
    end
end
end