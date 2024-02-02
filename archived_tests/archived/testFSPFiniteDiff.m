stoichMatrix = [1 0;-1 0;0 1;0 -1]';
fspTol = 1.0e-6;

kon = 1.0;
koff = 1.0;
kr = 1.0;
gamma = 0.1;

propStrings = cell(4, 1);
propStrings{1} = 'kon*(x1==0)';
propStrings{2} = 'koff*(x1==1)';
propStrings{3} = 'kr*(x1==1)';
propStrings{4} = 'gamma*x2';

par_names = {'kon', 'koff', 'kr', 'gamma'};

model = ssit.SrnModel(stoichMatrix, propStrings, par_names, []);

perturbationFactor = 1.0e-4;
n_max = 15;
x0 = [0 0]';
p0 = 1;
f_constr = @(x) [-x(1,:); -x(2,:); x(1,:); x(2,:)];
b_constr = [0 0 1 n_max]';
tspan = linspace(0, 10, 10);

parameters = {
    'kon' kon;
    'koff' koff;
    'kr' kr;
    'gamma' gamma
    };
out = ssit.sensitivity.adaptiveFspFiniteDiff(...
    model,...
    parameters,...
    perturbationFactor,...
    tspan,...
    x0,...
    p0,...
    fspTol,...
    f_constr, b_constr, 1);
pass = true;

