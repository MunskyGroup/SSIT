function [means, covs] = pgfMoments(sg, G)
% [means, covs] = ssit.pgf.pgfMoments(sg, G) means and covariances from the
% PGF node values G (N x Nt):
%   E[x_i]          = dG/dz_i (1),
%   E[x_i(x_j - d_ij)] = d2G/(dz_i dz_j) (1),
%   cov(i,j)        = d2G/(dz_i dz_j)(1) + d_ij E[x_i] - E[x_i] E[x_j].
% The derivatives are those of the sparse-grid interpolant at z = 1 (which
% equal the derivative-matrix rows when z = 1 is a grid node).
%   means : d x Nt,   covs : d x d x Nt.
d = sg.d;
Nt = size(G,2);
one = ones(1,d);
means = zeros(d, Nt);
for i = 1:d
    o = zeros(1,d); o(i) = 1;
    means(i,:) = ssit.pgf.interpolate(sg, G, one, o);
end
covs = zeros(d, d, Nt);
for i = 1:d
    for j = i:d
        o = zeros(1,d); o(i) = o(i) + 1; o(j) = o(j) + 1;
        f2 = ssit.pgf.interpolate(sg, G, one, o);
        v = f2 - means(i,:).*means(j,:);
        if i == j
            v = v + means(i,:);
        end
        covs(i,j,:) = reshape(v, 1, 1, Nt);
        covs(j,i,:) = covs(i,j,:);
    end
end
end
