function Phi = smolyakBasis(X, mi, domain, derivOrders)
% Phi = ssit.pgf.smolyakBasis(X, mi, domain, derivOrders) evaluates the
% Smolyak Chebyshev basis (or a partial derivative of it) at the points X.
%
%   X           : M x d evaluation points (e.g. the grid itself).
%   mi          : multi-index of the grid (ssit.pgf.smolyakMultiIndex).
%   domain      : 2 x d domain of the grid (ssit.pgf.checkDomain).
%   derivOrders : 1 x d derivative order per dimension, entries 0, 1 or 2
%                 (default zeros).  E.g. [0 1 0] = d/dz2, [1 1 0] =
%                 d2/dz1dz2, [2 0 0] = d2/dz1^2.
%
%   Phi         : M x N matrix, N = number of grid nodes.  With
%                 W = inv(smolyakBasis(grid, mi, domain)), Phi*W*G evaluates
%                 the requested derivative of the interpolant of the grid
%                 values G at X.
%
% One function replaces the interpolation, derivative and Hessian basis
% builders of Edward's utils.jl.  The column order is the same: blocks in
% multi-index order, and inside a block the LAST dimension's polynomial index
% varies fastest (rowwise_kron), so inverse/derivative matrices match the
% Julia ones entry by entry.
arguments
    X (:,:) double
    mi (:,:) double
    domain (2,:) double
    derivOrders (1,:) double = zeros(1, size(X,2))
end
[M, d] = size(X);
if size(mi,2) ~= d || numel(derivOrders) ~= d
    error('ssit:pgf:basisSize', 'X, the multi-index and derivOrders must all have d = %d columns.', d);
end
dom = ssit.pgf.checkDomain(domain, d);

% Map X to [-1,1] per dimension.
Xn = zeros(M, d);
for dim = 1:d
    if dom(1,dim) == dom(2,dim)
        Xn(:,dim) = 0;
    else
        Xn(:,dim) = 2*(X(:,dim) - dom(2,dim)) / (dom(1,dim) - dom(2,dim)) - 1;
    end
end

% Per level and dimension: the polynomial columns that are new at that level.
nLevels = max(mi(:));
newCols = cell(nLevels, d);
for lev = 1:nLevels
    m = ssit.pgf.smolyakLevelSize(lev);
    if lev == 1
        first = 1;
    else
        first = ssit.pgf.smolyakLevelSize(lev-1) + 1;
    end
    for dim = 1:d
        P = ssit.pgf.chebyshevBasis(m-1, Xn(:,dim), derivOrders(dim));
        newCols{lev,dim} = P(:, first:end);
    end
end

blocks = cell(1, size(mi,1));
for j = 1:size(mi,1)
    block = newCols{mi(j,1),1};
    for dim = 2:d
        block = rowKron(block, newCols{mi(j,dim),dim});
    end
    blocks{j} = block;
end
Phi = [blocks{:}];

% Chain rule for the map to [-1,1].
scale = prod((2 ./ (dom(1,:) - dom(2,:))) .^ derivOrders);
Phi = Phi * scale;
end

function C = rowKron(A, B)
% Row-wise Kronecker product with the column index of B varying fastest:
% C(:, (a-1)*q + b) = A(:,a) .* B(:,b).
[N, p] = size(A);
q = size(B,2);
C = reshape(reshape(B, N, q, 1) .* reshape(A, N, 1, p), N, p*q);
end
