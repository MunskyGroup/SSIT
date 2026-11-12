function G = fspToPgf(fspSolutions, z, speciesIdx)
% G = ssit.pgf.fspToPgf(fspSolutions, z) evaluates FSP solutions as PGFs,
% G(z) = sum_x p(x) prod_i z_i^x_i, at the points z (N x d).
% G = ssit.pgf.fspToPgf(fspSolutions, z, speciesIdx) PGF of the marginal of
% the species speciesIdx (one per column of z), e.g. the PGF species of a
% model with enumerated species (compare with sol.Gtot).
%   fspSolutions : SSIT FSP solutions (obj.Solutions.fsp), a cell array with
%                  one entry per time whose .p is an ssit.FspVector, or a
%                  single such entry.
%   G            : N x Nt.
% Used to compare the PGF solver with FSP.
if ~iscell(fspSolutions)
    fspSolutions = {fspSolutions};
end
if nargin < 3
    speciesIdx = 1:size(z,2);
end
Nt = numel(fspSolutions);
G = zeros(size(z,1), Nt);
for it = 1:Nt
    data = fspSolutions{it}.p.data;
    X = data.subs - 1;                 % states (subs are 1-based)
    p = data.vals;
    Zp = ones(size(z,1), size(X,1));
    for i = 1:size(z,2)
        Zp = Zp .* (z(:,i) .^ (X(:,speciesIdx(i))'));   % parentheses needed: ' and .^ share precedence
    end
    G(:,it) = Zp * p;
end
end
