function r = calibrateSSALogL(model, theta, levels, nReps, lossFun, kdeOpts, enforceIndependence)
% ssit.parest.calibrateSSALogL - Mean and std of the SSA log-likelihood
% estimate at each number of runs in levels (SSIT.ssaLikelihoodConvergence),
% with the KDE settings that the sequential ABC uses at that number of
% runs (ssit.parest.kdeArgsForRuns).
%
%   r = ssit.parest.calibrateSSALogL(model, theta, levels, nReps, lossFun, ...
%       kdeOpts, enforceIndependence)
%
% r has fields std, mean, values (numel(levels) x nReps), secondsPerEval.
r = struct('std', nan(1,numel(levels)), 'mean', nan(1,numel(levels)), ...
    'values', nan(numel(levels), nReps), 'secondsPerEval', nan(1,numel(levels)));
for k = 1:numel(levels)
    kdeArgs = ssit.parest.kdeArgsForRuns(kdeOpts, levels(k));
    res = model.ssaLikelihoodConvergence('nSimsList', levels(k), 'nReps', nReps, ...
        'lossFun', lossFun, 'pars', theta, 'enforceIndependence', enforceIndependence, ...
        'compareFSP', false, 'makePlot', false, kdeArgs{:});
    r.std(k) = res.std; r.mean(k) = res.mean;
    r.values(k,:) = res.values; r.secondsPerEval(k) = res.secondsPerEval;
end
end
