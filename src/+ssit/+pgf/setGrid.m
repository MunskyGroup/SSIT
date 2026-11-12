function prob = setGrid(prob, lower)
% prob = ssit.pgf.setGrid(prob, lower) builds the sparse grid on the domain
% [lower_j, 1] of each PGF species and assembles the PGF operators.
%   prob  : from ssit.pgf.prepare,
%   lower : 1 x dY lower ends in [-1, 1).
% Errors ssit:pgf:secondOrderDomain if a species with a quadratic
% propensity term does not have lower = -1.
model = prob.model;
Y = model.pgfIdx;
dY = numel(Y);
lower = lower(:)';
if numel(lower) ~= dY
    error('ssit:pgf:options', 'lowerBound must have one entry per PGF species (%d).', dY);
end
if any(lower < -1) || any(lower >= 1)
    error('ssit:pgf:options', 'pgfOptions.lowerBound must lie in [-1, 1).');
end
bad = model.secondOrder(Y) & lower > -1 + 1e-12;
if any(bad)
    error('ssit:pgf:secondOrderDomain', ...
        ['Species %s appear quadratically in a propensity (e.g. x*(x-1)). The PGF equation is then ', ...
         'second order and only well posed on the full domain [-1, 1] for these species: set their ', ...
         'pgfOptions.lowerBound to -1 (or use lowerBound = ''auto''), or enumerate them if their ', ...
         'counts are small.'], strjoin(model.species(Y(bad)), ', '));
end
prob.sg = ssit.pgf.buildGrid(prob.mu, [ones(1,dY); lower]);
prob.ops = ssit.pgf.assembleOperators(model, prob.sg);
prob.lower = lower;
end
