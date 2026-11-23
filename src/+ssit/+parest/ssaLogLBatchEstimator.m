function [logL, state] = ssaLogLBatchEstimator(model, theta, nTotal, state, lossFun, kdeOpts, enforceIndependence)
% ssit.parest.ssaLogLBatchEstimator - SSA-based log-likelihood estimate
% that can be refined by adding more SSA runs.
%
%   [logL, state] = ssit.parest.ssaLogLBatchEstimator(model, theta, nTotal, state, lossFun, kdeOpts)
%
% Runs enough new SSA trajectories at parameters theta to bring the pooled
% number of runs per time point up to nTotal, pools them with those already
% in state, and returns logL = -computeLossFunctionSSA(lossFun) for the
% pooled trajectories. Pass state = [] for a fresh estimate.
%
% Inputs:
%   model   -- SSIT model with solutionScheme 'SSA', loaded data, numeric
%              fittingOptions.modelVarsToFit, tSpan already including
%              initialTime, and fittingOptions.logPrior = [] (the prior is
%              handled by the caller).
%   theta   -- values of the parameters model.fittingOptions.modelVarsToFit
%              (linear space).
%   nTotal  -- SSA runs per time point after this call.
%   state   -- struct with fields npe (runs per time point so far) and sol
%              (pooled SSA solution), or [].
%   lossFun -- likelihood-type loss for computeLossFunctionSSA (e.g.,
%              'kde_loglik').
%   kdeOpts -- struct with kdeBandwidth, kdeMinBandwidth, kdeProbFloor,
%              kdeBiasCorrection, and optionally fixed widths
%              bandwidthTable / bandwidthLevels (see
%              ssit.parest.kdeArgsForRuns).
%   enforceIndependence -- (true) if true, each time point uses its own
%              block of nTotal trajectories (nTotal*numTimes runs in
%              all); if false, the same nTotal trajectories are used at
%              every time point (numTimes-fold fewer runs).
if nargin < 7, enforceIndependence = true; end
if isempty(state)
    state = struct('npe', 0, 'sol', []);
end
nAdd = nTotal - state.npe;
if nAdd > 0
    M = model;
    M.parameters(M.fittingOptions.modelVarsToFit,2) = num2cell(theta(:));
    M.ssaOptions.Nexp = 1;
    if enforceIndependence
        M.ssaOptions.nSimsPerExpt = nAdd;
    else
        M.ssaOptions.nSimsPerExpt = NaN;
        M.ssaOptions.Nsims = nAdd;
    end
    sol = M.solve(returnType='soln');
    state.sol = ssit.ssa.poolSSABatches(state.sol, state.npe, sol, nAdd);
    state.npe = nTotal;
end
M = model;
M.Solutions = state.sol;
M.ssaOptions.Nexp = 1;
M.ssaOptions.nSimsPerExpt = state.npe;
kdeArgs = ssit.parest.kdeArgsForRuns(kdeOpts, state.npe);
logL = -M.computeLossFunctionSSA(lossFun, [], enforceIndependence, true, kdeArgs{:});
end
