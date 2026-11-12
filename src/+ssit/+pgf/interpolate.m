function V = interpolate(sg, G, X, derivOrders)
% V = ssit.pgf.interpolate(sg, G, X) evaluates the sparse-grid interpolant
% of the node values G (N x k, one column per function, e.g. per time) at
% the points X (M x d).  Returns M x k.
%
% V = ssit.pgf.interpolate(sg, G, X, derivOrders) evaluates a partial
% derivative instead (derivOrders: 1 x d, entries 0/1/2), e.g. the means
% E[x_i] = dG/dz_i at X = ones(1,d) when the grid has no z = 1 node.
arguments
    sg (1,1) struct
    G (:,:) double
    X (:,:) double
    derivOrders (1,:) double = zeros(1, sg.d)
end
if size(G,1) ~= sg.N
    error('ssit:pgf:interpSize', 'G must have one row per grid node (%d).', sg.N);
end
Phi = ssit.pgf.smolyakBasis(X, sg.multiIndex, sg.domain, derivOrders);
V = Phi * (sg.invInterp * G);
end
