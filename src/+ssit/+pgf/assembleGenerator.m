function L = assembleGenerator(ops, t, theta)
% L = ssit.pgf.assembleGenerator(ops, t, theta) PGF generator at time t and
% parameters theta, so that dG/dt = L G with G = [G_1; ...; G_nE] (one block
% of N grid values per enumerated state).
%   Single block (no enumerated species): dense N x N.
%   Several blocks: sparse (nE N) x (nE N) with dense N x N blocks.
N = ops.N;
nE = ops.nE;
C = ssit.pgf.termCoefficients(ops, t, theta);
if nE == 1
    L = reshape(ops.Mflat * C, N, N);
    return
end

model = ops.model;
blocks = cell(nE, nE);
for q = 1:ops.nTerms
    target = model.rxn(ops.termReaction(q)).target;
    for e = 1:nE
        c = C(q,e);
        if c == 0 || target(e) == 0
            continue
        end
        f = target(e);
        blocks{f,e} = addTo(blocks{f,e}, c * ops.Mgain(:,:,q));
        blocks{e,e} = addTo(blocks{e,e}, -c * ops.Mloss(:,:,q));
    end
end

[ii, jj] = ndgrid(1:N, 1:N);
[fi, ei] = find(~cellfun(@isempty, blocks));
nb = numel(fi);
I = zeros(N*N, nb); J = zeros(N*N, nb); V = zeros(N*N, nb);
for b = 1:nb
    I(:,b) = (fi(b)-1)*N + ii(:);
    J(:,b) = (ei(b)-1)*N + jj(:);
    V(:,b) = blocks{fi(b), ei(b)}(:);
end
L = sparse(I(:), J(:), V(:), nE*N, nE*N);
end

function B = addTo(B, M)
if isempty(B)
    B = M;
else
    B = B + M;
end
end
