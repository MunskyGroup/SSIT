function app = runMoments(app)
% this function calculates the mean, variance and covariance from moment
% closure equations

inputs = app.ReactionsTabOutputs.inputs; % Pulls the input equations specified
pars = app.ReactionsTabOutputs.parameters; % Pulls the parameters specified
stoich = app.ReactionsTabOutputs.stoichMatrix; % Pulls the stoichiometry from the reactions
propens = app.ReactionsTabOutputs.propensities; % Pulls the propensity functions

propensSym=str2sym(propens);
propensParam=subs(propensSym,pars(:,1),pars(:,2));
inputSym=str2sym(inputs(:,2));
PropensSymFinal=subs(propensParam,inputs(:,1),inputSym);

syms x1 x2 x3 real
w=PropensSymFinal'; % propensity matrix in symbolic form
S=stoich; % stoich matrix
x = [x1;x2;x3]; % species
n = size(x,1); % number of species

vvt=[];
try
    ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS');
    
    x0 = eval(app.MomentsInitCondField.Value);   % Pulls the inital conditions specified
    v0 = [x0(1) x0(2) x0(3),...
        x0(1)+x0(1)^2 x0(1)*x0(2) x0(1)*x0(3),...
        x0(2)+x0(2)^2 x0(2)*x0(3),...
        x0(3)+x0(3)^2]'; % create initial conditions
    momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
    tspan = eval(app.MomentsPrintTimesField.Value); % Pulls the array of times
    [t,vvt] = ode45(momOde,tspan,v0); % solve the ode from time 0 to 10
catch
    warndlg('model propensity functions are invalid for moment calculations','error')
end

app.MomentTabOutputs.moments = vvt;
end