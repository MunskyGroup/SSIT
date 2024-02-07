stoichMatrix = [1]';
fspTol = 1.0e-4;

wstr = ssit.SrnModel.processPropensityStrings({'1'});

n_max = 10;
x0 = 0;
p0 = 1;
f_constr = @(x) [-x(:); x(:)];
b_constr = [0 n_max]';

propensities = cell(1,1);
for i = 1:1
    propensities{i} = ssit.Propensity.createFromString(wstr{i}, stoichMatrix(:,i), i);
end

tspan = linspace(0, 10, 10);

[Outputs, b_constr] = ssit.fsp.adaptiveFspSolve(tspan, x0, p0, stoichMatrix, ...
                        propensities, fspTol, f_constr, b_constr, 0);

pass = true;
for i = 1:5
    err = 0;
    for j = 1:nnz(Outputs{i}.p.data)
        err = err + abs(Outputs{i}.p.data.values(j) - poisspdf(j-1, 1.0*tspan(i)));
    end    
    if (err > fspTol)
        pass = false;
        break;
    end
end

assert(pass);
