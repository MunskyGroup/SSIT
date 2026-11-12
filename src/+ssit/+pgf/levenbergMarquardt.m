function [x, fval, info] = levenbergMarquardt(fp, x0, fitOptions)
% [x, fval, info] = ssit.pgf.levenbergMarquardt(fp, x0, fitOptions)
% minimises -logL = fp.scale * ||r(exp(x))||^2 over the log-parameters x
% (natural log, order fp.fitIdx) for the fit set up by ssit.pgf.fitSetup,
% with r = ssit.pgf.residuals.  Port of Edward's Levenberg-Marquardt
% (PGFfromEdward/test_smolyak_grdusp1_lm_steady.jl): forward-difference
% Jacobian in log space, Marquardt damping lambda*diag(J'J), step cap, box
% around the start, lambda / 5 after an accepted step and * 5 after a
% rejected one.  No toolbox needed.
% Called by SSIT.maximizeLikelihood(..., 'levenbergMarquardt').
%
% fitOptions (struct, optional; optimset fields are honoured):
%   MaxIter    (100)   iterations (one Jacobian each),
%   TolFun     (1e-5)  stop when the relative decrease of -logL is below,
%   TolX       (1e-7)  stop when the step norm is below,
%   Display    ('iter') 'iter', 'final' or 'off',
%   lambda0    (1e-2)  initial damping,
%   stepCap    (2)     largest change of any log-parameter per step,
%   logBound   (log(1e4)) box |x - x0| <= logBound,
%   fdStep     (1e-3)  forward-difference step in log space,
%   maxRejected (40)   rejected trial steps per iteration before stopping,
%   UseParallel (false) parfor over the Jacobian columns,
%   computeFIM (true)  Jacobian at the solution for info.fimLog.
%
% Returns x (log-parameters), fval = -logL at x, and info with fields
%   iterations, nfev, converged, message, trace ([iteration nfev seconds
%   -logL x']), J, fimLog (Gauss-Newton Fisher information of the
%   log-parameters, 2*scale*J'*J; meaningful for loss 'gls').
arguments
    fp (1,1) struct
    x0 (:,1) double
    fitOptions = struct()
end
o = struct('MaxIter', 100, 'TolFun', 1e-5, 'TolX', 1e-7, 'Display', 'iter', ...
    'lambda0', 1e-2, 'stepCap', 2, 'logBound', log(1e4), 'fdStep', 1e-3, ...
    'maxRejected', 40, 'UseParallel', false, 'computeFIM', true);
names = fieldnames(o);
for i = 1:numel(names)
    if isstruct(fitOptions) && isfield(fitOptions, names{i}) && ~isempty(fitOptions.(names{i}))
        o.(names{i}) = fitOptions.(names{i});
    end
end
verbose = strcmpi(o.Display, 'iter');
resFun = @(x) ssit.pgf.residuals(fp, exp(x));
tStart = tic;
x = x0;
r = resFun(x);
f = sum(r.^2);
nfev = 1;
if ~isfinite(f)
    error('ssit:pgf:fit', 'The PGF solve fails at the initial parameters.');
end
lambda = o.lambda0;
trace = [0, nfev, 0, fp.scale*f, x.'];
converged = false;
message = 'maximum number of iterations reached';
J = [];
it = 0;
if verbose
    fprintf('  iter        -logL      lambda    max|step|   nfev   time(s)\n');
    fprintf('  %4d  %11.6g\n', 0, fp.scale*f);
end
for it = 1:o.MaxIter
    J = jacobian(resFun, x, r, o.fdStep, o.UseParallel);
    nfev = nfev + numel(x);
    if any(~isfinite(J(:)))
        message = 'PGF solve failed while computing the Jacobian';
        break
    end
    A = J.' * J;
    g = J.' * r;
    accepted = false;
    for attempt = 1:o.maxRejected
        step = -(A + lambda * diag(diag(A) + 1e-12)) \ g;
        mstep = max(abs(step));
        if mstep > o.stepCap
            step = step * (o.stepCap / mstep);
        end
        step = min(max(x + step, x0 - o.logBound), x0 + o.logBound) - x;
        rNew = resFun(x + step);
        nfev = nfev + 1;
        fNew = sum(rNew.^2);
        if fNew < f
            x = x + step;
            r = rNew;
            lambda = max(lambda / 5, 1e-12);
            accepted = true;
            trace(end+1, :) = [it, nfev, toc(tStart), fp.scale*fNew, x.']; %#ok<AGROW>
            if verbose
                fprintf('  %4d  %11.6g  %10.2e  %10.3g  %5d  %8.1f\n', it, fp.scale*fNew, lambda, ...
                    max(abs(step)), nfev, toc(tStart));
            end
            done = (f - fNew) < o.TolFun * f || norm(step) < o.TolX;
            f = fNew;
            break
        else
            lambda = lambda * 5;
        end
    end
    if ~accepted
        message = 'no decrease found (local minimum or step too small)';
        converged = true;
        break
    end
    if done
        message = 'relative decrease below TolFun or step below TolX';
        converged = true;
        break
    end
end
fval = fp.scale * f;
if o.computeFIM
    J = jacobian(resFun, x, r, o.fdStep, o.UseParallel);
    nfev = nfev + numel(x);
end
info = struct('iterations', it, 'nfev', nfev, 'converged', converged, 'message', message, ...
    'trace', trace, 'J', J, 'fimLog', 2 * fp.scale * (J.' * J), 'time', toc(tStart));
if ~strcmpi(o.Display, 'off')
    fprintf('Levenberg-Marquardt: %s after %d iterations (%d PGF solves, %.1f s), -logL = %.8g\n', ...
        message, it, nfev, info.time, fval);
end
end

function J = jacobian(resFun, x, r0, h, useParallel)
% Forward differences in log space.
P = numel(x);
J = zeros(numel(r0), P);
nWorkers = 0;
if useParallel
    nWorkers = Inf;
end
parfor (k = 1:P, nWorkers)
    xk = x;
    xk(k) = xk(k) + h;
    J(:,k) = (resFun(xk) - r0) / h;
end
end
