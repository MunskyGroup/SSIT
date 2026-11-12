function sg = buildGrid(mu, domain, opts)
% sg = ssit.pgf.buildGrid(mu, domain) builds a Smolyak sparse grid and the
% matrices needed to apply derivatives of the PGF at its nodes.
%
%   mu     : scalar (isotropic) or 1 x d vector (anisotropic) Smolyak layer.
%   domain : 2 x d box, row 1 = upper bounds (normally 1), row 2 = lower
%            bounds (either orientation accepted).
%   Name-value options:
%     hessian (false) : also build all second-derivative matrices H.
%
% Returns struct sg with fields
%   d, mu, domain (normalised), z (N x d nodes), multiIndex, N,
%   invInterp  : N x N, maps node values G to Chebyshev coefficients,
%   D          : N x N x d, D(:,:,i)*G = dG/dz_i at the nodes,
%   H          : N x N x d x d (empty unless requested),
%                H(:,:,a,b)*G = d2G/(dz_a dz_b) at the nodes,
%   oneIndex   : row of z equal to all ones (0 if absent, which happens
%                when max(mu) < d; use ssit.pgf.interpolate at z = 1 then).
%
% Port of the grid / derivative-matrix construction in Edward's Julia code
% (PGFfromEdward/utils.jl); reference data: tests/test_data/pgf.
arguments
    mu (1,:) double {mustBeInteger, mustBeNonnegative}
    domain (2,:) double
    opts.hessian (1,1) logical = false
end
d = size(domain,2);
dom = ssit.pgf.checkDomain(domain, d);
[z, mi] = ssit.pgf.smolyakGrid(mu, dom);
N = size(z,1);

Phi0 = ssit.pgf.smolyakBasis(z, mi, dom);
invInterp = inv(Phi0);

D = zeros(N, N, d);
for i = 1:d
    orders = zeros(1,d); orders(i) = 1;
    D(:,:,i) = ssit.pgf.smolyakBasis(z, mi, dom, orders) * invInterp;
end

H = [];
if opts.hessian
    H = zeros(N, N, d, d);
    for a = 1:d
        for b = a:d
            orders = zeros(1,d);
            orders(a) = orders(a) + 1;
            orders(b) = orders(b) + 1;
            H(:,:,a,b) = ssit.pgf.smolyakBasis(z, mi, dom, orders) * invInterp;
            H(:,:,b,a) = H(:,:,a,b);
        end
    end
end

oneIndex = find(all(abs(z - 1) <= 1e-12, 2), 1);
if isempty(oneIndex)
    oneIndex = 0;
end

sg = struct('d', d, 'mu', mu, 'domain', dom, 'z', z, 'multiIndex', mi, ...
    'N', N, 'invInterp', invInterp, 'D', D, 'H', H, 'oneIndex', oneIndex);
end
