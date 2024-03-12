import ssit.*

stoichMatrix = [1 0;-1 0;0 1;0 -1]';
fspTol = 1.0e-4;

kon = 1.0;
koff = 1.0;
kr = 1.0;
gamma = 0.1;

propensities = cell(4, 1);
propensityDerivatives = cell(4, 4);

propensities{1} = Propensity.createAsTimeInvariant(@(x) kon*x(1,:)==0, stoichMatrix(:,1), 1);
propensities{2} = Propensity.createAsTimeInvariant(@(x) koff*x(1,:)==1, stoichMatrix(:,2), 2);
propensities{3} = Propensity.createAsTimeInvariant(@(x) kr*x(1,:)==1, stoichMatrix(:,3), 3);
propensities{4} = Propensity.createAsTimeInvariant(@(x) gamma*x(2,:), stoichMatrix(:,4), 4);

propensityDerivatives{1,1} = Propensity.createAsTimeInvariant(@(x) x(1,:)==0, stoichMatrix(:,1), 1);
propensityDerivatives{2,2} = Propensity.createAsTimeInvariant(@(x) x(1,:)==1, stoichMatrix(:,2), 2);
propensityDerivatives{3,3} = Propensity.createAsTimeInvariant(@(x) x(1,:)==1, stoichMatrix(:,3), 3);
propensityDerivatives{4,4} = Propensity.createAsTimeInvariant(@(x) x(2,:), stoichMatrix(:,4), 4);

n_max = 15;
x0 = [0 0]';
p0 = 1;
s0 = [0 0 0 0]';

f_constr = @(x) [-x(1,:); -x(2,:); x(1,:); x(2,:)];
b_constr = [0 0 1 n_max]';

tspan = linspace(0, 10, 10);

[Outputs, b_constraints] = ssit.sensitivity.adaptiveFspForwardSens(tspan, x0,...
    p0, s0, stoichMatrix, propensities, propensityDerivatives, fspTol, f_constr, b_constr,...
    0, 1, 1.0E-4, 1.0E-14);

pass = true;
assert(pass);
