function G = propagate(prob, theta, checkStability)
% G = ssit.pgf.propagate(prob, theta) block PGF at the times prob.tSpan for
% the parameter values theta (default prob.theta), on the grid set by
% ssit.pgf.setGrid.  Returns G (nE*N x Nt) = [G_1; ...; G_nE].
% G = ssit.pgf.propagate(prob, theta, checkStability) overrides
% prob.opts.checkStability (the eigenvalue check costs one dense eig).
%
% The initial PGF is the steady state of L(tSpan(1)) if
% prob.init.steadyState, else built from initialCondition / initialProbs.
% Errors: ssit:pgf:unstable, ssit:pgf:steadyState, ssit:pgf:integration.
% (No arguments block: this function has a nested function.)
if nargin < 2
    theta = [];
end
if nargin < 3
    checkStability = [];
end
if isempty(prob.sg)
    error('ssit:pgf:noGrid', 'Call ssit.pgf.setGrid (or ssit.pgf.chooseGrid) first.');
end
if isempty(theta)
    theta = prob.theta;
end
if isempty(checkStability)
    checkStability = prob.opts.checkStability;
end
model = prob.model;
sg = prob.sg;
ops = prob.ops;
tSpan = prob.tSpan;
opts = prob.opts;
cache = struct('C', [], 'L', []);
Lfun = @generatorAt;

if checkStability
    checkTimes = tSpan(1);
    if model.timeDependent
        checkTimes = unique([tSpan(1), tSpan(end)]);
    end
    for t = checkTimes
        L = Lfun(t);
        maxRe = max(real(eig(full(L))));
        if maxRe > 1e-8 * max(1, norm(L, 1))
            involved = model.secondOrder | model.mixedOrder;
            error('ssit:pgf:unstable', ...
                ['The discretised PGF operator is unstable (max real eigenvalue %.3g at t = %g). ', ...
                 'This happens for bimolecular reactions on a truncated z-domain. Widen the domain ', ...
                 '(pgfOptions.lowerBound toward -1) for species %s, lower pgfOptions.mu, or enumerate ', ...
                 'species with small counts.'], maxRe, t, strjoin(model.species(involved), ', '));
        end
    end
end

if prob.init.steadyState
    G0 = ssit.pgf.steadyState(ops, sg, tSpan(1), theta);
else
    G0 = blockInitialPgf(model, sg, prob.init.x0, prob.init.p0);
end
G = integrate(Lfun, model.timeDependent, tSpan, G0, opts);

    function Lt = generatorAt(tt)
        % L(tt), reassembled only when the term coefficients change (e.g.
        % once per segment for a step input): the coefficients are cheap,
        % the assembly is not.
        if model.timeDependent || isempty(cache.L)
            Ct = ssit.pgf.termCoefficients(ops, tt, theta);
            if ~isequal(Ct, cache.C)
                cache.C = Ct;
                cache.L = ssit.pgf.assembleGenerator(ops, tt, theta);
            end
        end
        Lt = cache.L;
    end
end

% -------------------------------------------------------------------------
function G0 = blockInitialPgf(model, sg, x0, p0)
% Initial block PGF: G_e(z) = sum_m p0(m) 1{x0_E(:,m) = e} z^x0_Y(:,m).
N = sg.N;
nE = size(model.enumStates, 2);
M = size(x0, 2);
if numel(p0) ~= M
    if M == 1
        p0 = 1;
    else
        error('ssit:pgf:initialCondition', ...
            'initialProbs must have one entry per column of initialCondition (%d).', M);
    end
end
G0 = zeros(N*nE, 1);
for m = 1:M
    if isempty(model.enumIdx)
        e = 1;
    else
        [~, e] = ismember(x0(model.enumIdx, m)', model.enumStates', 'rows');
    end
    rows = (e-1)*N + (1:N);
    G0(rows) = G0(rows) + ssit.pgf.initialPgf(sg.z, x0(model.pgfIdx, m), p0(m));
end
end

function G = integrate(Lfun, timeDependent, tSpan, G0, opts)
% Stiff implicit integration of dG/dt = L(t) G, restarted at timeBreaks.
% Within each segment (a, b] the generator is evaluated right-continuously
% at a (L(a + delta)), so a switch such as (t > a) acts from the start of the
% segment, while a steady state computed at t = a uses L(a).
Nt = numel(tSpan);
G = zeros(numel(G0), Nt);
G(:,1) = G0;
if Nt == 1
    return
end
solver = str2func(opts.odeSolver);
breaks = opts.timeBreaks(:)';
edges = unique([tSpan(1), breaks(breaks > tSpan(1) & breaks < tSpan(end)), tSpan(end)]);
g = G0;
for e = 1:numel(edges)-1
    a = edges(e); b = edges(e+1);
    odeOpts = odeset('RelTol', opts.relTol, 'AbsTol', opts.absTol);
    if timeDependent
        delta = 1e-10 * max(1, b - a);
        Lseg = @(t) Lfun(max(t, a + delta));
        rhs = @(t,y) Lseg(t) * y;
        odeOpts = odeset(odeOpts, 'Jacobian', @(t,y) Lseg(t));
    else
        L = Lfun(a);
        rhs = @(t,y) L * y;
        odeOpts = odeset(odeOpts, 'Jacobian', L);
    end
    % |G(z)| <= 1 on the domain (|z| <= 1), so |G| > 2 means the solution
    % has gone unstable (e.g. at a bad trial point during a fit): stop
    % early instead of grinding through tiny steps.
    odeOpts = odeset(odeOpts, 'OutputFcn', @stopIfUnbounded);
    % A failed step is reported as an error below (the fitting routines
    % treat it as an infinite objective), so silence the solver's warning.
    warnState = warning('off', ['MATLAB:', opts.odeSolver, ':IntegrationTolNotMet']);
    cleanup = onCleanup(@() warning(warnState));
    odeSol = solver(rhs, [a b], g, odeOpts);
    clear cleanup
    if odeSol.x(end) < b - 1e-12 * max(1, abs(b))
        error('ssit:pgf:integration', ...
            'The PGF integration stopped at t = %g before reaching t = %g.', odeSol.x(end), b);
    end
    inSeg = find(tSpan > a & tSpan <= b);
    if ~isempty(inSeg)
        G(:,inSeg) = deval(odeSol, tSpan(inSeg));
    end
    g = deval(odeSol, b);
end
if ~all(isfinite(G(:)))
    error('ssit:pgf:integration', 'The PGF integration produced non-finite values.');
end
end

function stop = stopIfUnbounded(~, y, flag)
stop = false;
if isempty(flag) && ~isempty(y)
    stop = ~all(isfinite(y(:))) || max(abs(y(:))) > 2;
end
end
