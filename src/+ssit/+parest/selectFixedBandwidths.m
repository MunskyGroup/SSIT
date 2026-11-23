function bandwidthTable = selectFixedBandwidths(model, theta, levels, nReps, lossFun, kdeOpts, enforceIndependence)
% ssit.parest.selectFixedBandwidths - Choose fixed KDE kernel widths for
% each number of SSA runs, fitted time and observed species.
%
%   bandwidthTable = ssit.parest.selectFixedBandwidths(model, theta, levels, nReps, ...
%       lossFun, kdeOpts, enforceIndependence)
%
% For each N = levels(k), runs nReps independent SSA simulations at the
% parameters theta, lets computeLossFunctionSSA choose the widths with
% kdeOpts.kdeBandwidth (normally 'lcv'), and stores the median over the
% repeats: bandwidthTable{k} is a numTimes x d matrix (d = number of
% observed species) that can be passed as kdeBandwidth. Using the same
% widths for every parameter set removes the (costly) bandwidth search
% from each likelihood evaluation, and the shared smoothing makes the
% log-likelihood differences between parameter sets less noisy.
bandwidthTable = cell(1, numel(levels));
for k = 1:numel(levels)
    M = model;
    M.ssaOptions.Nexp = 1;
    if enforceIndependence
        M.ssaOptions.nSimsPerExpt = levels(k);
    else
        M.ssaOptions.nSimsPerExpt = NaN;
        M.ssaOptions.Nsims = levels(k);
    end
    H = [];
    for r = 1:nReps
        [~, ~, h] = M.computeLossFunctionSSA(lossFun, theta, enforceIndependence, false, ...
            kdeBandwidth=kdeOpts.kdeBandwidth, kdeMinBandwidth=kdeOpts.kdeMinBandwidth, ...
            kdeProbFloor=kdeOpts.kdeProbFloor, kdeBiasCorrection=kdeOpts.kdeBiasCorrection);
        H = cat(3, H, h);
    end
    bandwidthTable{k} = median(H, 3);
end
end
