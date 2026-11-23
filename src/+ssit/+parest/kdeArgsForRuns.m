function args = kdeArgsForRuns(kdeOpts, n)
% ssit.parest.kdeArgsForRuns - KDE name-value arguments for
% SSIT.computeLossFunctionSSA (and ssaLikelihoodConvergence) when the SSA
% likelihood is estimated from n runs per time point.
%
%   args = ssit.parest.kdeArgsForRuns(kdeOpts, n)
%
% kdeOpts fields: kdeBandwidth, kdeMinBandwidth, kdeProbFloor,
% kdeBiasCorrection, optionally kdeOutlierFraction / kdeOutlierSupport
% (outlier-robust likelihood), and optionally bandwidthMode ('lcv', 'local',
% 'fixed'), bandwidthTable (cell, one numTimes x d matrix of reference
% widths per entry of bandwidthLevels) and localFactors. The reference
% widths used are those of the largest level <= n (the first level if n
% is smaller than all levels).
%   'lcv'   (or no table): kdeBandwidth as given (full search).
%   'local': kdeBandwidth = 'lcv' searching localFactors x reference.
%   'fixed': kdeBandwidth = reference widths.
args = {'kdeMinBandwidth', kdeOpts.kdeMinBandwidth, 'kdeProbFloor', kdeOpts.kdeProbFloor, ...
    'kdeBiasCorrection', kdeOpts.kdeBiasCorrection};
if isfield(kdeOpts, 'kdeOutlierFraction')
    args = [args, {'kdeOutlierFraction', kdeOpts.kdeOutlierFraction}];
end
if isfield(kdeOpts, 'kdeOutlierSupport')
    args = [args, {'kdeOutlierSupport', kdeOpts.kdeOutlierSupport}];
end
mode = 'lcv';
if isfield(kdeOpts, 'bandwidthMode'), mode = kdeOpts.bandwidthMode; end
hasTable = isfield(kdeOpts, 'bandwidthTable') && ~isempty(kdeOpts.bandwidthTable);
if strcmp(mode, 'lcv') || ~hasTable
    args = [args, {'kdeBandwidth', kdeOpts.kdeBandwidth}];
    return
end
k = find(kdeOpts.bandwidthLevels <= n, 1, 'last');
if isempty(k), k = 1; end
ref = kdeOpts.bandwidthTable{k};
switch mode
    case 'fixed'
        args = [args, {'kdeBandwidth', ref}];
    case 'local'
        args = [args, {'kdeBandwidth', 'lcv', 'kdeLcvCenter', ref, ...
            'kdeLcvFactors', kdeOpts.localFactors}];
end
end
