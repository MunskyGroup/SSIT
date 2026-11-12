function [r, H] = residuals(fp, pars)
% r = ssit.pgf.residuals(fp, pars) weighted residual vector of the PGF fit
% set up by ssit.pgf.fitSetup, at the fitted parameters pars (linear scale,
% order of fp.fitIdx; [] = current values), such that
%   logL = -fp.scale * sum(r.^2)
% (see ssit.pgf.logLikelihood).  A log prior adds the entries
% sqrt(-logPrior/scale), which assumes each entry of logPrior(pars) is <= 0
% (e.g. Gaussian or log-normal priors).  If the PGF solve fails, every
% entry is Inf.
% [r, H] = ... also returns the model values compared with the data
% ((M*nG) x numel(fp.fitTimes), see ssit.pgf.comparisonValues).
theta = fp.theta0;
if nargin > 1 && ~isempty(pars)
    theta(fp.fitIdx) = pars(:);
else
    pars = theta(fp.fitIdx);
end
H = [];
r = inf(fp.nRes, 1);
if any(~isfinite(theta)) || any(theta < 0)
    return
end
prob = fp.prob;
try
    G = ssit.pgf.propagate(prob, theta, false);
catch err
    if startsWith(err.identifier, 'ssit:pgf:')
        return
    end
    rethrow(err)
end
nE = size(prob.model.enumStates, 2);
H = ssit.pgf.comparisonValues(G, prob.sg.N, nE, fp.tCol, fp.Agg, fp.Phi);
nk = size(H, 1);
for k = 1:numel(fp.tCol)
    dk = H(:,k) - fp.data(k).m;
    W = fp.data(k).W;
    if isvector(W) && ~isscalar(W)
        rk = W .* dk;
    else
        rk = W * dk;
    end
    r((k-1)*nk + (1:nk)) = rk;
end
if ~isempty(fp.logPrior)
    lp = fp.logPrior(pars(:));
    r(numel(fp.tCol)*nk + 1:end) = sqrt(max(0, -lp(:) / fp.scale));
end
if any(~isfinite(r))
    r(:) = inf;
end
end
