function logL = logLikelihood(fp, pars)
% logL = ssit.pgf.logLikelihood(fp, pars) PGF (approximate) log-likelihood
% of the data, plus the log prior, for the fit set up by ssit.pgf.fitSetup:
%   logL = -fp.scale * sum(r.^2),  r = ssit.pgf.residuals(fp, pars).
% pars are the fitted parameters on the linear scale (order fp.fitIdx; []
% = current values).  Returns -Inf if the PGF solve fails.
%   loss 'gls' / 'weighted' : Gaussian approximation of the log-likelihood
%                             of the empirical PGF (up to a constant).
%   loss 'mse'              : minus the mean squared PGF residual.
if nargin < 2
    pars = [];
end
r = ssit.pgf.residuals(fp, pars);
logL = -fp.scale * sum(r.^2);
if ~isfinite(logL)
    logL = -Inf;
end
end
