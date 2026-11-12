function [z, mi] = smolyakGrid(mu, domain)
% [z, mi] = ssit.pgf.smolyakGrid(mu, domain) builds the Smolyak sparse grid
% on Chebyshev-Gauss-Lobatto nodes, reproducing
% SmolyakApprox.smolyak_grid(chebyshev_gauss_lobatto, d, mu, domain).
%
%   mu     : scalar (isotropic) or 1 x d vector (anisotropic) layer, see
%            ssit.pgf.smolyakMultiIndex.
%   domain : 2 x d, row 1 = upper bounds, row 2 = lower bounds (either
%            orientation is accepted, see ssit.pgf.checkDomain).
%
%   z      : N x d grid nodes, one node per row.
%   mi     : multi-index (one row per tensor block).
%
% Level i of the 1-D rule has m(i) = 1 (i = 1) or 2^(i-1)+1 nodes; each
% level contributes only the nodes that are new relative to level i-1.  The
% nodes of each block are the tensor product of these new nodes, with the
% first dimension varying fastest.
arguments
    mu (1,:) double {mustBeInteger, mustBeNonnegative}
    domain (2,:) double
end
d = size(domain,2);
dom = ssit.pgf.checkDomain(domain, d);
mi = ssit.pgf.smolyakMultiIndex(d, mu);

nLevels = max(mi(:));
newNodes = cell(nLevels,1);
previous = [];
for lev = 1:nLevels
    nodes = ssit.pgf.chebyshevNodes(ssit.pgf.smolyakLevelSize(lev));
    newNodes{lev} = setdiff(nodes, previous, 'stable');
    previous = nodes;
end

blocks = cell(size(mi,1),1);
for j = 1:size(mi,1)
    block = newNodes{mi(j,1)};
    for dim = 2:d
        next = newNodes{mi(j,dim)};
        n1 = size(block,1);
        block = [repmat(block, numel(next), 1), kron(next(:), ones(n1,1))];
    end
    blocks{j} = block;
end
x = vertcat(blocks{:});

% Scale from [-1,1] to the domain.
z = dom(2,:) + (1 + x) .* (dom(1,:) - dom(2,:)) * 0.5;
end
