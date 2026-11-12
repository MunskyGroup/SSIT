function [sol, cache] = solve(src, opts, cache)
% sol = ssit.pgf.solve(src) solves the chemical master equation of an SSIT
% model through its probability generating function (PGF)
%   G(z,t) = E[prod_j z_j^y_j(t)]
% on a Smolyak sparse grid in z (Edward's method, PGFfromEdward/).  Small,
% bounded species (gene states) are enumerated: one PGF block per reachable
% enumerated state (see ssit.pgf.translateModel and pgfOptions.enumeratedSpecies).
%
%   src  : SSIT object (or struct with the same fields): species,
%          propensityFunctions, stoichiometry, parameters, inputExpressions,
%          initialCondition, initialProbs, tSpan, (fspOptions.initApproxSS).
%   opts : options (default src.pgfOptions, else ssit.pgf.defaultOptions);
%          missing fields take their default values.
%
% Returns struct sol with fields
%   T_array    : 1 x Nt output times (src.tSpan),
%   G          : nE*N x Nt block PGF [G_1; ...; G_nE] at the grid nodes,
%   Gtot       : N x Nt PGF of the PGF species (sum over enumerated states),
%   means      : d x Nt,  covs : d x d x Nt   (all species, model order),
%   enumeratedSpecies, pgfSpecies : names,
%   enumeratedStates : dE x nE,  enumeratedProbabilities : nE x Nt,
%   grid       : sparse-grid struct over the PGF species (ssit.pgf.buildGrid),
%   lowerBound : 1 x dY lower ends of the z-domain used,
%   species, parameters, options, solveTime.
%
% [sol, cache] = ssit.pgf.solve(src, opts, cache) reuses the setup of a
% previous solve: the translated model (ssit.pgf.translateModel), the sparse
% grid and its operators, and (lowerBound = 'auto') the z-domain.  They are
% reused when the model structure is unchanged -- species, propensities,
% stoichiometry, parameter names, inputs, initial condition and the
% enumeration options -- and tSpan stays within the time range of the
% cached solve; parameter values, tolerances and the initial distribution
% may differ.  The grid is also kept when mu and lowerBound (numeric) are
% unchanged; with lowerBound = 'auto' the cached domain is kept as long as
% it still suits the new solution (see ssit.pgf.chooseGrid).  Pass
% cache = [] to rebuild everything.
%
% Use SSIT.solve(solver='pgf') to call this from an SSIT model (the cache
% is kept in the hidden SSIT property pgfCache).
% Requirements: propensities polynomial (degree <= 2) in the PGF species;
% species with a quadratic term need lowerBound = -1 (set by 'auto').
arguments
    src
    opts = []
    cache = []
end
tic;
opts = ssit.pgf.mergeOptions(src, opts);
tSpan = src.tSpan(:)';
key = cacheKey(src, opts);
reused = '';
if ~isempty(cache) && isequal(cache.key, key) && ...
        tSpan(1) >= cache.tRange(1) && tSpan(end) <= cache.tRange(2)
    prob = ssit.pgf.prepare(src, opts, tSpan, cache.model);
    reused = 'model';
    if isequal(prob.mu, cache.mu)
        prob.sg = cache.sg;
        prob.ops = cache.ops;
        prob.lower = cache.lower;
        reused = 'model and grid';
    end
else
    prob = ssit.pgf.prepare(src, opts, tSpan);
    cache = struct('key', {key}, 'tRange', tSpan([1 end]));
end
[prob, G] = ssit.pgf.chooseGrid(prob);
if strcmp(reused, 'model and grid') && ~isequal(prob.lower, cache.lower)
    reused = 'model';                 % the domain no longer fitted: grid rebuilt
end
cache.model = prob.model;
cache.mu = prob.mu;
cache.sg = prob.sg;
cache.ops = prob.ops;
cache.lower = prob.lower;
model = prob.model;
[means, covs, P, Gtot] = ssit.pgf.computeMoments(model, prob.sg, G);
sol = struct('T_array', prob.tSpan, 'G', G, 'Gtot', Gtot, 'means', means, 'covs', covs, ...
    'enumeratedProbabilities', P, 'grid', prob.sg, 'lowerBound', prob.lower);
sol.species = model.species;
sol.enumeratedSpecies = model.species(model.enumIdx);
sol.pgfSpecies = model.species(model.pgfIdx);
sol.enumeratedStates = model.enumStates;
sol.parameters = model.parameterValues';
sol.options = prob.opts;
sol.solveTime = toc;
if prob.opts.verbose
    fprintf('PGF solve: enumerated {%s} (%d states), N = %d nodes, lowerBound = [%s], %.2f s\n', ...
        strjoin(sol.enumeratedSpecies, ', '), size(model.enumStates,2), sol.grid.N, ...
        num2str(sol.lowerBound, '%.4g '), sol.solveTime);
    if ~isempty(reused)
        fprintf('           (reused the %s of the previous solve)\n', reused);
    end
end
end

function key = cacheKey(src, opts)
% Everything the translated model depends on, apart from parameter values
% and the probe times (see ssit.pgf.translateModel).
inputs = {};
if (isstruct(src) && isfield(src, 'inputExpressions')) || (isobject(src) && isprop(src, 'inputExpressions'))
    inputs = src.inputExpressions;
end
key = {src.species(:), cellstr(string(src.propensityFunctions(:))), src.stoichiometry, ...
    src.parameters(:,1), inputs, src.initialCondition, ...
    opts.enumeratedSpecies, opts.maxEnumeratedCount, opts.maxEnumeratedStates};
end
