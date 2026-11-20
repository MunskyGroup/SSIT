function [smpl, accept, value, info] = sequentialMHSample(start, nsamples, estimator, levels, sigma, opts)
% ssit.parest.sequentialMHSample - Metropolis-Hastings with a noisy,
% refinable log-likelihood estimate and a sequential accept/reject test.
%
%   [smpl, accept, value, info] = ssit.parest.sequentialMHSample( ...
%       start, nsamples, estimator, levels, sigma, proprnd=@(x)..., ...)
%
% At each step, with current state x and proposal y:
%   1. Draw u ~ U(0,1) and set the threshold
%        tau = log(u) - [logprior(y) - logprior(x)].
%      Exact MH accepts y exactly when LLD = logL(y) - logL(x) >= tau.
%   2. Estimate logL at x and at y from levels(1) runs each (fresh,
%      independent runs for both, so that their biases at equal N cancel),
%      giving LLDhat with standard deviation sqrt(2)*sigma(1).
%   3. If LLDhat - tau >  z*sqrt(2)*sigma(k), accept;
%      if LLDhat - tau < -z*sqrt(2)*sigma(k), reject.
%   4. Otherwise add runs to both estimates (pooled with the earlier runs)
%      to reach levels(k+1) and repeat. At the last level, the plain MH
%      rule is used: accept if LLDhat >= tau.
%
% Inputs:
%   start     -- (1 x p) initial state (row vector).
%   nsamples  -- number of samples to return.
%   estimator -- @(x, n, state) returning [logL, state]: estimate of logL
%                at x refined to n runs; state = [] starts a fresh
%                estimate (see ssit.parest.ssaLogLBatchEstimator).
%   levels    -- increasing numbers of runs, e.g. [300 600 ... 9600].
%   sigma     -- std of a single logL estimate at each level (calibrated).
% Name-value options:
%   proprnd           -- (required) proposal sampler @(x) (symmetric).
%   logprior          -- (@(x)0) log prior density of the state.
%   z                 -- (2.5) confidence multiplier.
%   thin, burnin      -- (1, 0) as in metropolisHastingsSample.
%   recalibrateEvery  -- (Inf) recalibrate sigma after every K-th
%                        acceptance, using recalibrateFun.
%   recalibrateFun    -- ([]) @(x) returning a new sigma vector at x.
%   saveFileName      -- ('') if set, progress is saved after every step
%                        and an interrupted run resumes from this file.
%   verbose           -- (true) print one line per step.
%
% Outputs:
%   smpl   -- (nsamples x p) chain.
%   accept -- acceptance rate over all steps after burn-in.
%   value  -- (nsamples x 1) estimated logL + logprior of each sample
%             (estimated with the number of runs used at that step, so
%             values are noisy and not all at the same N).
%   info   -- struct: decisionLevel (index into levels of the level at
%             which each step was decided; 0 = rejected without
%             simulation, e.g. zero prior), decidedByMHRule (true where
%             the last level was reached without a confident decision),
%             accepted, simsPerParameterSet, sigmaHistory (one row per
%             calibration), recalibratedAtStep, bestfound, bestValue.
arguments
    start (1,:) double
    nsamples (1,1) double {mustBePositive, mustBeInteger}
    estimator function_handle
    levels (1,:) double {mustBePositive}
    sigma (1,:) double {mustBeNonnegative}
    opts.proprnd function_handle
    opts.logprior function_handle = @(x) 0
    opts.z (1,1) double {mustBePositive} = 2.5
    opts.thin (1,1) double {mustBePositive, mustBeInteger} = 1
    opts.burnin (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    opts.recalibrateEvery (1,1) double {mustBePositive} = Inf
    opts.recalibrateFun = []
    opts.saveFileName = ''
    opts.verbose (1,1) logical = true
end
if ~isfield(opts, 'proprnd')
    error('SSIT:sequentialMHSample:noProposal', 'A proposal sampler proprnd is required.');
end
if numel(sigma) ~= numel(levels)
    error('SSIT:sequentialMHSample:sigma', 'sigma must have one entry per level.');
end
if any(diff(levels) <= 0)
    error('SSIT:sequentialMHSample:levels', 'levels must be increasing.');
end

nSteps = opts.burnin + nsamples*opts.thin;
p = numel(start);
saveFile = char(opts.saveFileName);
if ~isempty(saveFile) && isfile(saveFile)
    S = load(saveFile);
    smpl = S.smpl; value = S.value; info = S.info; x = S.x; lpx = S.lpx;
    sigma = S.sigma; iStart = S.iStep + 1; nAcc = S.nAcc; nAccAll = S.nAccAll;
    fprintf('Resuming sequential MH from step %d of %d (%s).\n', iStart, nSteps, saveFile);
else
    smpl = nan(nsamples, p);
    value = nan(nsamples, 1);
    info = struct('decisionLevel', zeros(nSteps,1), 'decidedByMHRule', false(nSteps,1), ...
        'accepted', false(nSteps,1), 'simsPerParameterSet', zeros(nSteps,1), ...
        'sigmaHistory', sigma(:)', 'recalibratedAtStep', 0, ...
        'bestfound', start, 'bestValue', -Inf, 'levels', levels, 'z', opts.z);
    x = start;
    lpx = opts.logprior(x);
    if ~isfinite(lpx)
        error('SSIT:sequentialMHSample:badStart', 'The log prior is not finite at the start.');
    end
    iStart = 1; nAcc = 0; nAccAll = 0;
end
curValue = NaN;

for iStep = iStart:nSteps
    y = opts.proprnd(x);
    lpy = opts.logprior(y);
    acc = false; k = 0; logLx = NaN; logLy = NaN;
    if isfinite(lpy)
        tau = log(rand) - (lpy - lpx);
        stX = []; stY = [];
        for k = 1:numel(levels)
            [logLx, stX] = estimator(x, levels(k), stX);
            [logLy, stY] = estimator(y, levels(k), stY);
            if ~isfinite(logLy)
                acc = false; break
            end
            margin = (logLy - logLx) - tau;
            s = opts.z*sqrt(2)*sigma(k);
            if margin > s
                acc = true; break
            elseif margin < -s
                acc = false; break
            elseif k == numel(levels)
                acc = margin >= 0;            % plain MH rule at the last level
                info.decidedByMHRule(iStep) = true;
            end
        end
        info.simsPerParameterSet(iStep) = levels(k);
    end
    info.decisionLevel(iStep) = k;
    info.accepted(iStep) = acc;

    if acc
        x = y; lpx = lpy; curValue = logLy + lpy;
        nAccAll = nAccAll + 1;
        if iStep > opts.burnin, nAcc = nAcc + 1; end
        if curValue > info.bestValue
            info.bestValue = curValue; info.bestfound = x;
        end
        if isfinite(opts.recalibrateEvery) && mod(nAccAll, opts.recalibrateEvery) == 0 ...
                && ~isempty(opts.recalibrateFun)
            sigma = opts.recalibrateFun(x);
            info.sigmaHistory(end+1,:) = sigma(:)';
            info.recalibratedAtStep(end+1) = iStep;
        end
    elseif isfinite(logLx)
        curValue = logLx + lpx;
    end

    if iStep > opts.burnin && mod(iStep - opts.burnin, opts.thin) == 0
        j = (iStep - opts.burnin)/opts.thin;
        smpl(j,:) = x;
        value(j) = curValue;
    end
    if opts.verbose
        if k == 0
            fprintf('step %d/%d: rejected (zero prior)\n', iStep, nSteps);
        else
            fprintf('step %d/%d: %s at N = %d%s; LLDhat = %.3g\n', iStep, nSteps, ...
                ternary(acc,'accepted','rejected'), levels(k), ...
                ternary(info.decidedByMHRule(iStep),' (MH rule)',''), logLy - logLx);
        end
    end
    if ~isempty(saveFile)
        save(saveFile, 'smpl', 'value', 'info', 'x', 'lpx', 'sigma', 'iStep', 'nAcc', 'nAccAll');
    end
end
accept = nAcc / max(1, nSteps - opts.burnin);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
