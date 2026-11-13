function [p, info] = validatedQSS(jacR, jacB, v0)
% ssit.fsp.validatedQSS - Validated quasi-steady-state (QSS) distribution
% of a truncated FSP generator.
%
%   [p, info] = ssit.fsp.validatedQSS(jacR, jacB)
%   [p, info] = ssit.fsp.validatedQSS(jacR, jacB, v0)
%
% Inputs:
%   * jacR - (n x n sparse) FSP generator restricted to the retained states
%            (columns are "from" states).
%   * jacB - (m x n sparse) rates from each retained state into the m
%            FSP sinks/constraints.
%   * v0   - (optional, n x 1) starting vector for eigs, e.g., a previously
%            validated QSS padded with zeros. Ignored if empty or wrong size.
%
% Outputs:
%   * p    - (n x 1) nonnegative QSS distribution, sum(p) = 1.
%   * info - struct with fields
%       ok       true if p passed validation
%       method   method that produced p
%       lambda   QSS eigenvalue (0 for the reflecting solve)
%       escape   total rate of probability leaving the FSP from p
%                (= -lambda for the Perron vector); use as the FSP error
%                indicator for the steady state
%       residual ||jacR*p - lambda*p||_1
%       message  reason for the last rejection (empty if first try worked)
%
% Method:
%   The QSS is the Perron vector of jacR: its eigenvalue lambda <= 0 is real
%   and rightmost, and because every other eigenvalue has a smaller real
%   part it is also the eigenvalue of smallest magnitude. The QSS is
%   therefore computed by shift-invert eigs ('smallestabs'), which converges
%   in a few iterations even when the spectrum of jacR is many orders of
%   magnitude wider than |lambda| (the case where eigs(...,'largestreal')
%   fails to converge and returns a meaningless vector). Each candidate is
%   accepted only if eigs converged, the vector is finite and single-signed,
%   the eigen-residual is small, and lambda = -escape (the columns of a
%   generator sum to zero, so this holds for any true Perron vector).
%   Candidates are tried in the order:
%     1. eigs(jacR,1,'smallestabs')
%     2. eigs(jacR,1,sigma) with a tiny negative shift sigma
%     3. eigs(jacR,1,'largestreal')
%     4. Reflecting linear solve (last resort): probability leaving each
%        state through the FSP boundary is returned to that state, giving a
%        proper generator Q (columns sum to zero); the stationary vector of
%        Q is obtained from one sparse direct solve of Q*p = 0, sum(p) = 1.
%        This approximates the same steady state and converges to it as the
%        FSP is expanded; escape(p) is again the error indicator.
%   Warnings from these attempts are suppressed only for their duration and
%   the caller's warning state is restored afterwards.
arguments
    jacR
    jacB
    v0 = []
end

n = size(jacR,1);
escapeRates = full(sum(jacB,1)).';
info = struct('ok',false,'method','','lambda',NaN,'escape',NaN, ...
    'residual',NaN,'message','');

% Trivial case.
if n == 1
    p = 1;
    info.ok = true; info.method = 'trivial';
    info.lambda = full(jacR); info.escape = escapeRates; info.residual = 0;
    return
end

% Suppress (only) the warnings that the attempts below may raise, and
% restore the caller's warning state on exit.
warnState = warning;
restoreWarnings = onCleanup(@()warning(warnState));
warnIds = {'MATLAB:eigs:NotAllEigsConverged','MATLAB:eigs:SigmaNearExactEig', ...
    'MATLAB:eigs:IllConditionedA','MATLAB:nearlySingularMatrix', ...
    'MATLAB:singularMatrix','MATLAB:illConditionedMatrix'};
for k = 1:numel(warnIds)
    warning('off',warnIds{k});
end

normA = normest1(jacR);

eigsOpts = struct();
if ~isempty(v0) && numel(v0) == n && any(v0(:)~=0) && all(isfinite(v0(:)))
    eigsOpts.v0 = full(v0(:));
end

% 1-3: validated Perron vector.
targets = {'smallestabs', -1e-12*normA, 'largestreal'};
targetNames = {'smallestabs','shift','largestreal'};
for k = 1:numel(targets)
    try
        [v,lam,flag] = eigs(jacR,1,targets{k},eigsOpts);
    catch err
        info.message = sprintf('%s: eigs error (%s)', targetNames{k}, err.message);
        continue
    end
    [p, cand] = checkPerron(v, lam, flag, jacR, escapeRates, normA);
    if cand.ok
        info = cand;
        info.method = targetNames{k};
        return
    end
    info.message = sprintf('%s: %s', targetNames{k}, cand.message);
end

% 4: reflecting linear solve (last resort).
try
    Q = jacR + spdiags(escapeRates, 0, n, n);   % columns now sum to zero
    [~,kRow] = max(abs(diag(Q)));              % replace one equation by sum(p)=1
    Q(kRow,:) = 1;
    b = zeros(n,1); b(kRow) = 1;
    p = Q\b;
    if all(isfinite(p)) && min(p) >= -1e-8*max(abs(p))
        p = max(real(p),0); p = p/sum(p);
        info.ok = true;
        info.method = 'reflect';
        info.lambda = 0;
        info.escape = escapeRates.'*p;
        info.residual = norm(Q*p - b, 1);
        return
    end
    info.message = [info.message, '; reflect: non-finite or negative solution'];
catch err
    info.message = [info.message, '; reflect: ', err.message];
end

% Nothing worked: uniform distribution, flagged as not ok.
p = ones(n,1)/n;
info.method = 'uniform';
info.escape = escapeRates.'*p;
end

function [p, info] = checkPerron(v, lam, flag, jacR, escapeRates, normA)
% Validate a candidate Perron vector of jacR.
info = struct('ok',false,'method','','lambda',lam,'escape',NaN, ...
    'residual',NaN,'message','');
p = [];
if flag ~= 0
    info.message = 'eigs did not converge'; return
end
if ~isfinite(lam) || any(~isfinite(v))
    info.message = 'non-finite eigenvalue or eigenvector'; return
end
lam = real(lam);
v = real(v);
v = v*sign(sum(v));                         % fix the arbitrary overall sign
if ~(max(v) > 0) || min(v) < -1e-8*max(abs(v))
    info.message = 'eigenvector is not single-signed'; return
end
p = max(v,0); p = p/sum(p);
info.lambda = lam;
info.residual = norm(jacR*p - lam*p, 1);
info.escape = escapeRates.'*p;
roundoff = 1e3*eps*normA;
if info.residual > 1e-8*normA + roundoff
    info.message = sprintf('large eigen-residual (%.2e)', info.residual); return
end
if abs(info.escape + lam) > 1e-6*max(abs(lam),info.escape) + roundoff
    info.message = sprintf('eigenvalue (%.2e) inconsistent with escape rate (%.2e)', ...
        lam, info.escape); return
end
info.ok = true;
end
