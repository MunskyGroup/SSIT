function mi = smolyakMultiIndex(d, mu)
% mi = ssit.pgf.smolyakMultiIndex(d, mu) returns the Smolyak multi-index
% (one row per tensor block, entries are 1-based levels) for a d-dimensional
% sparse grid, with the same row order as SmolyakApprox.generate_multi_index:
%
%   scalar mu, d > 1 : isotropic grid of layer mu, built recursively
%                      (rows [i, base_j] with i varying fastest).
%   vector mu (1 x d): anisotropic grid; candidates are enumerated with the
%                      first index varying fastest and kept when
%                      sum(row) <= d + max(mu) and row <= mu + 1.
%
% A scalar mu gives the same index (same rows, same order) as the vector
% repmat(mu,1,d); this was checked for d <= 4, mu <= 3 and is asserted in
% tests/pgfGridTest.m.  Both rules are kept to mirror SmolyakApprox.jl.
arguments
    d (1,1) double {mustBeInteger, mustBePositive}
    mu (1,:) double {mustBeInteger, mustBeNonnegative}
end
if isscalar(mu)
    mi = isotropicIndex(d, mu);
else
    if numel(mu) ~= d
        error('ssit:pgf:muSize', 'mu has %d entries but the grid has %d dimensions.', numel(mu), d);
    end
    mi = anisotropicIndex(d, mu);
end
end

function mi = isotropicIndex(d, mu)
if d == 1
    mi = (1:mu+1)';
    return
end
base = isotropicIndex(d-1, mu);
mi = zeros(nchoosek(mu+d, d), d);
pos = 0;
for j = 1:size(base,1)
    rowSum = sum(base(j,:));
    for i = 1:mu+1
        if rowSum + i <= d + mu
            pos = pos + 1;
            mi(pos,:) = [i, base(j,:)];
        end
    end
end
mi = mi(1:pos,:);
end

function mi = anisotropicIndex(d, mu)
maxMu = max(mu);
B = maxMu + 1;
k = (0:B^d-1)';               % zero-based linear index, first subscript fastest
subs = zeros(numel(k), d);
r = k;
for dim = 1:d
    subs(:,dim) = mod(r, B) + 1;
    r = floor(r / B);
end
keep = sum(subs,2) <= d + maxMu & all(subs <= mu + 1, 2);
mi = subs(keep,:);
end
