function ops = assembleOperators(model, sg)
% ops = ssit.pgf.assembleOperators(model, sg) builds the constant matrices
% of the (block) PGF generator on the grid sg (over the PGF species).  For
% each term q = (reaction r, falling-factorial order k):
%
%   Mgain_q = diag(z^(s_Y + k)) * d^k   (into block e + s_E),
%   Mloss_q = diag(z^k) * d^k           (out of block e),
%
% so that, with C = ssit.pgf.termCoefficients(ops, t, theta) (nTerms x nE),
%
%   L(target_r(e), e) += C(q,e) * Mgain_q,    L(e, e) -= C(q,e) * Mloss_q.
%
% Use ssit.pgf.assembleGenerator to form L.  Derivative matrices d^k are
% exact derivatives of the sparse-grid interpolant (ssit.pgf.smolyakBasis),
% including mixed and second derivatives.
%
% Returns struct ops with fields
%   Mgain, Mloss : N x N x nTerms,
%   Mflat        : N^2 x nTerms, columns (Mgain - Mloss)(:)  (single-block shortcut),
%   termReaction, termRow, N, nE, nTerms, model.
arguments
    model (1,1) struct
    sg (1,1) struct
end
N = sg.N;
z = sg.z;
nTerms = sum(arrayfun(@(x) size(x.k,1), model.rxn));
Mgain = zeros(N, N, nTerms);
Mloss = zeros(N, N, nTerms);
termReaction = zeros(nTerms,1);
termRow = zeros(nTerms,1);
derivCache = containers.Map('KeyType','char','ValueType','any');

q = 0;
for r = 1:model.R
    sY = model.S(model.pgfIdx, r)';
    for row = 1:size(model.rxn(r).k,1)
        k = model.rxn(r).k(row,:);
        key = sprintf('%d_', k);
        if ~isKey(derivCache, key)
            if all(k == 0)
                derivCache(key) = eye(N);
            else
                derivCache(key) = ssit.pgf.smolyakBasis(z, sg.multiIndex, sg.domain, k) * sg.invInterp;
            end
        end
        Dk = derivCache(key);
        q = q + 1;
        Mgain(:,:,q) = prod(z.^(sY + k), 2) .* Dk;
        Mloss(:,:,q) = prod(z.^k, 2) .* Dk;
        termReaction(q) = r;
        termRow(q) = row;
    end
end

nE = size(model.enumStates, 2);
ops = struct('Mgain', Mgain, 'Mloss', Mloss, ...
    'Mflat', reshape(Mgain - Mloss, N*N, nTerms), ...
    'N', N, 'nE', nE, 'nTerms', nTerms, ...
    'termReaction', termReaction, 'termRow', termRow, 'model', model);
end
