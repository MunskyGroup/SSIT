function G = steadyState(ops, sg, t, theta)
% G = ssit.pgf.steadyState(ops, sg, t, theta) stationary PGF of the
% generator frozen at time t: L(t) G = 0 with total probability one,
%   sum_e G_e(z = 1) = 1,
% where G_e(1) is evaluated through the sparse-grid interpolant (so a z = 1
% grid node is not required).  Solved as the overdetermined system
% [L; normalisation] G = [0; 1] by QR.  For a single block on a grid with a
% z = 1 node this reproduces Edward's row-replacement solve.
% Errors ssit:pgf:steadyState if the result is not finite or L G is not ~0.
N = ops.N;
nE = ops.nE;
L = ssit.pgf.assembleGenerator(ops, t, theta);
phi = ssit.pgf.smolyakBasis(ones(1, sg.d), sg.multiIndex, sg.domain) * sg.invInterp;   % 1 x N
normRow = repmat(phi, 1, nE);
A = [full(L); normRow];
b = [zeros(nE*N, 1); 1];
% Economy QR without column pivoting: A can be very ill conditioned (the
% operators are non-normal, cond ~ 1e13 for GR-DUSP1 at all-ones), and the
% rank decision of A\b would then return a "basic" solution with zeros.
[Q, R] = qr(A, 0);
warnState = warning('off', 'MATLAB:nearlySingularMatrix');
cleanup = onCleanup(@() warning(warnState));
G = R \ (Q' * b);
clear cleanup
res = max(abs(L * G));
if ~all(isfinite(G)) || res > 1e-8 * max(1, norm(L, 1)) * max(1, max(abs(G)))
    error('ssit:pgf:steadyState', ...
        ['The steady-state PGF could not be found (max |L G| = %.3g). The model may have no unique ', ...
         'steady state at t = %g (e.g. a conserved total among the PGF species).'], res, t);
end
end
