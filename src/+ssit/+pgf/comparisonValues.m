function H = comparisonValues(G, N, nE, cols, Agg, Phi)
% H = ssit.pgf.comparisonValues(G, N, nE, cols, Agg, Phi) model PGF values
% that are compared with the data (see ssit.pgf.fitSetup):
%   H(:,k) = vec( Phi * G_blocks(:,:,cols(k)) * Agg' ),
% i.e. for each time column cols(k) of the block PGF G (nE*N x Nt) the
% blocks are summed within each group of enumerated states (Agg: nG x nE)
% and evaluated at the comparison points (Phi: M x N interpolation matrix,
% [] when the points are the grid nodes).  H is (M*nG) x numel(cols),
% group-major (all points of group 1, then group 2, ...).
Gb = reshape(G, N, nE, []);
nG = size(Agg, 1);
if isempty(Phi)
    M = N;
else
    M = size(Phi, 1);
end
H = zeros(M*nG, numel(cols));
for k = 1:numel(cols)
    Hk = Gb(:,:,cols(k)) * Agg.';
    if ~isempty(Phi)
        Hk = Phi * Hk;
    end
    H(:,k) = Hk(:);
end
end
