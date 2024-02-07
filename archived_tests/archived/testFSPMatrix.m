import ssit.Propensity

stoichMatrix = [1;-1;1]';

wstr = cell(3,1);
wstr{1} = '1.0';
wstr{2} = '1+cos(t)';
wstr{3} = '10.0/(1.0 + (1.0 + sin(t))*x1)';

nMax = 100;
x0 = 0;
f_constr = @(x) [-x(1,:); x(1,:)];
b_constr = [0 nMax]';

propensities = cell(3,1);
for i = 1:3
    propensities{i} = ssit.Propensity.createFromSym(str2sym(wstr{i}), stoichMatrix(:,i), i);
end
state_set = ssit.FiniteStateSet(x0, stoichMatrix);
state_set = state_set.expand(f_constr, b_constr);

FSPOperator = ssit.FspMatrix(propensities, state_set, length(b_constr));
%%
V = normrnd(0, 1, nMax+1+length(b_constr), 30);
pass = true;
for t = 0:0.1:1
    W1 = FSPOperator.multiply(t, V);
    At = GenerateAt(t, nMax);
    W2 = At*V;
    rel_err = norm(W1 - W2, 'fro')/norm(W1, 'fro');
    if (rel_err > eps)
        pass = false;
        break;
    end
end
assert(pass);
%%
function A = GenerateAt(t, n_max)

A1 = spdiags([-ones(n_max +3, 1) ones(n_max+3, 1)], [0 -1], n_max+3, n_max+3);

A2 = (1+cos(t))*spdiags([-ones(n_max+3, 1) ones(n_max+3,1)], [0 1], n_max+3, n_max+3);

v3 = 10.0./(1 + (1 + sin(t))*(0:n_max+2)');
A3 = spdiags([-v3 v3], [0 -1], n_max+3, n_max+3);
A = A1 + A2 + A3;
A(n_max+3, :) = A(n_max+3, :) + A(n_max+2, :);
A(n_max+2, :) = 0;
A(n_max+2, 1) = (1+cos(t));
A(:, n_max+2) = 0;
A(:, n_max+3) = 0;
end
