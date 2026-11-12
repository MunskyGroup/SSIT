function C = termCoefficients(ops, t, theta)
% C = ssit.pgf.termCoefficients(ops, t, theta) returns the nTerms x nE
% coefficients c_q(e, t, theta) of the PGF generator terms (see
% ssit.pgf.assembleOperators): for each reaction the propensity is
% evaluated on its small lattice of PGF-species values, for every
% enumerated state, and mapped to falling-factorial coefficients.
model = ops.model;
nE = ops.nE;
C = zeros(ops.nTerms, nE);
q0 = 0;
for r = 1:model.R
    x = model.rxn(r);
    nq = size(x.k,1);
    if nq == 0
        continue
    end
    v = x.w(x.Xfull, t, theta);
    v = real(v) + zeros(1, size(x.Xfull,2));   % complex zeros (checked in translateModel)
    C(q0+1:q0+nq, :) = x.A * reshape(v, size(x.X,2), nE);
    q0 = q0 + nq;
end
end
