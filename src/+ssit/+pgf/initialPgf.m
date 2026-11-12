function G0 = initialPgf(z, x0, p0)
% G0 = ssit.pgf.initialPgf(z, x0, p0) PGF of the initial distribution at
% the grid nodes z (N x d): G0(z) = sum_m p0(m) prod_i z_i^x0(i,m).
%   x0 : d x M initial states (SSIT initialCondition),
%   p0 : M x 1 probabilities (SSIT initialProbs; scalar 1 for a single state).
arguments
    z (:,:) double
    x0 (:,:) double
    p0 (:,1) double = 1
end
M = size(x0,2);
if numel(p0) ~= M
    if M == 1
        p0 = 1;
    else
        error('ssit:pgf:initialCondition', ...
            'initialProbs must have one entry per column of initialCondition (%d).', M);
    end
end
G0 = zeros(size(z,1),1);
for m = 1:M
    G0 = G0 + p0(m) * prod(z .^ (x0(:,m)'), 2);   % parentheses needed: ' and .^ share precedence
end
end
