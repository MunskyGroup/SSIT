function [logL, state] = ssaLogLBatchEstimator(model, theta, nTotal, state, lossFun, kdeOpts)
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
%              kdeBiasCorrection.
if isempty(state)
    state = struct('npe', 0, 'sol', []);
end
nAdd = nTotal - state.npe;
if nAdd > 0
    M = model;
    M.parameters(M.fittingOptions.modelVarsToFit,2) = num2cell(theta(:));
    M.ssaOptions.Nexp = 1;
    M.ssaOptions.nSimsPerExpt = nAdd;
    sol = M.solve(returnType='soln');
    state.sol = ssit.ssa.poolSSABatches(state.sol, state.npe, sol, nAdd);
    state.npe = nTotal;
end
M = model;
M.Solutions = state.sol;
M.ssaOptions.Nexp = 1;
M.ssaOptions.nSimsPerExpt = state.npe;
logL = -M.computeLossFunctionSSA(lossFun, [], true, true, ...
    kdeBandwidth=kdeOpts.kdeBandwidth, kdeMinBandwidth=kdeOpts.kdeMinBandwidth, ...
    kdeProbFloor=kdeOpts.kdeProbFloor, kdeBiasCorrection=kdeOpts.kdeBiasCorrection);
end
