function opts = defaultOptions()
% opts = ssit.pgf.defaultOptions() returns the default settings of the PGF
% solver (SSIT property pgfOptions, used by SSIT.solve(solver='pgf')).
%
% Grid
%   mu          Smolyak layer: scalar, or a vector with one entry per PGF
%               species (= species that are not enumerated) or per species.
%   lowerBound  Lower end of the z-domain [lower, 1] per PGF species: 'auto',
%               a scalar, or a vector (one entry per PGF species or per
%               species; values in [-1, 1)).
%               'auto' solves once on a provisional domain, then sets
%               lower_i = max(-1, 1 - autoScale / max_t E[x_i]) and
%               re-solves (autoPasses times).  Species whose propensities
%               are quadratic in them always get lower = -1 (required for
%               stability, see tests/test_data/pgf/README.md).
%   autoScale   Target decay of the PGF across the domain: for a Poisson
%               species G(lower) ~ exp(-autoScale).
%   autoPasses  Number of re-solves used by lowerBound = 'auto'.
%
% Time integration
%   odeSolver   'ode15s' (default) or 'ode23s'.  Use a stiff implicit
%               solver: the PGF operators can be strongly non-normal, so
%               expm/explicit methods are unreliable.
%   relTol, absTol  integrator tolerances.
%   timeBreaks  Times at which propensities jump (e.g. a step input); the
%               integration is restarted there.
%
% Enumerated (gene-type) species
%   enumeratedSpecies  'auto' (default), {} (none), or a cell array of
%               species names.  Enumerated species are not PGF variables:
%               their reachable states are listed explicitly and the PGF of
%               the other species is carried for each of them (blocks
%               G_e(z)).  Propensities may then be ANY function of the
%               enumerated species, and only need to be polynomial in the
%               others.  'auto' selects species bounded by
%                 (a) a nonnegative conservation law with a small total
%                     (e.g. offGene + onGene = 1), or
%                 (b) their own propensities (every reaction that increases
%                     the species vanishes above a small count, e.g.
%                     kon*(1-geneState)).
%   maxEnumeratedCount   largest count of an automatically enumerated species.
%   maxEnumeratedStates  cap on the number of enumerated states (error above).
%
% Initial condition
%   initialSteadyState  [] (follow fspOptions.initApproxSS), true or false:
%               start from the steady state of the generator at the
%               initial time instead of initialCondition.
%
% Checks
%   checkStability  Verify max real part of eig(L) <= 0 (at the first and
%               last time) before integrating; errors with advice otherwise.
%   verbose     Print the chosen domain and timings.
%
% Fitting to data (SSIT.maximizeLikelihood / computeLikelihood with
% solutionScheme = 'pgf'; see ssit.pgf.fitSetup)
%   loss        How the model PGF is compared with the empirical PGF of the
%               data, (1/n) sum_cells z^x, at the grid nodes:
%               'gls' (default)  Gaussian approximation of the likelihood
%                   of the empirical PGF (central limit theorem): residuals
%                   weighted by the inverse of the full covariance
%                   Cov(z_j^x, z_k^x)/n + errorFloor^2 * I estimated from
%                   the data.  Closest to the maximum-likelihood fit; logL
%                   is on the log-likelihood scale (usable for MCMC).
%               'weighted'  same, but only the variances (diagonal).
%               'mse'  unweighted mean squared PGF residual (Edward's
%                   objective); logL = -MSE.
%   errorFloor  Standard deviation added to the sampling noise of each
%               PGF value to account for the error of the PGF solution
%               itself ('gls' and 'weighted').  'auto' estimates it once,
%               at the initial parameters, from the change of the solution
%               when every mu is increased by one (at least relTol).
%   With lowerBound = 'auto', the domain of observed species is set from the
%   data means (and that of unobserved PGF species from the model at the
%   initial parameters), then kept fixed during the fit.
opts = struct( ...
    'enumeratedSpecies', 'auto', ...
    'maxEnumeratedCount', 10, ...
    'maxEnumeratedStates', 200, ...
    'initialSteadyState', [], ...
    'mu', 3, ...
    'lowerBound', 'auto', ...
    'autoScale', 2, ...
    'autoPasses', 2, ...
    'odeSolver', 'ode15s', ...
    'relTol', 1e-6, ...
    'absTol', 1e-9, ...
    'timeBreaks', [], ...
    'checkStability', true, ...
    'verbose', false, ...
    'loss', 'gls', ...
    'errorFloor', 'auto');
end
