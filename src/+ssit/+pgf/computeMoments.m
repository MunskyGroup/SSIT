function [means, covs, P, Gtot] = computeMoments(model, sg, G)
% [means, covs, P, Gtot] = ssit.pgf.computeMoments(model, sg, G) moments of
% ALL species (original order) from the block PGF G = [G_1; ...; G_nE]
% (nE*N x Nt):
%   P(e,t)  = G_e(1)                      probability of enumerated state e,
%   Gtot    = sum_e G_e                   PGF of the PGF species (N x Nt),
%   E[y], cov(y)  from Gtot               (ssit.pgf.pgfMoments),
%   E[x_E]  = sum_e e P(e),   cov(x_E) = sum_e e e' P(e) - E E',
%   E[x_E,i y_j] = sum_e e_i dG_e/dz_j(1).
%   means : d x Nt,  covs : d x d x Nt,  P : nE x Nt.
N = sg.N;
nE = size(model.enumStates, 2);
Nt = size(G, 2);
d = model.d;
Y = model.pgfIdx;
Eidx = model.enumIdx;
dY = numel(Y);
Gb = reshape(G, N, nE, Nt);
Gtot = reshape(sum(Gb, 2), N, Nt);
[mY, cY] = ssit.pgf.pgfMoments(sg, Gtot);

means = zeros(d, Nt);
covs = zeros(d, d, Nt);
means(Y,:) = mY;
covs(Y,Y,:) = cY;

one = ones(1, dY);
P = zeros(nE, Nt);
dG = zeros(dY, nE, Nt);            % dG_e/dz_j at z = 1
for e = 1:nE
    Ge = reshape(Gb(:,e,:), N, Nt);
    P(e,:) = ssit.pgf.interpolate(sg, Ge, one);
    if ~isempty(Eidx)
        for j = 1:dY
            o = zeros(1, dY); o(j) = 1;
            dG(j,e,:) = reshape(ssit.pgf.interpolate(sg, Ge, one, o), 1, 1, Nt);
        end
    end
end
if isempty(Eidx)
    return
end

X = model.enumStates;              % dE x nE
mE = X * P;
means(Eidx,:) = mE;
for it = 1:Nt
    covs(Eidx,Eidx,it) = X * (P(:,it) .* X') - mE(:,it) * mE(:,it)';
    EXY = X * dG(:,:,it)';         % dE x dY
    covs(Eidx,Y,it) = EXY - mE(:,it) * mY(:,it)';
    covs(Y,Eidx,it) = covs(Eidx,Y,it)';
end
end
