function figs = plotSSAStateFit(model, state, lossFun, kdeOpts, value, figurePrefix, enforceIndependence)
% ssit.parest.plotSSAStateFit - Plot data vs model marginals and means/stds
% vs time for the pooled SSA runs in an estimator state (see
% ssit.parest.ssaLogLBatchEstimator), replacing earlier plots of this kind.
% Used by SSIT.runSequentialABC(plotNewBest=true) when the chain finds a
% new best parameter set.
%
%   figs = ssit.parest.plotSSAStateFit(model, state, lossFun, kdeOpts, value, figurePrefix)
%
% value is the estimated log posterior (shown in the titles); if
% figurePrefix is not empty, the figures are saved as
% <figurePrefix>_marginals.(fig/png) and <figurePrefix>_meanStd.(fig/png).
if nargin < 7, enforceIndependence = true; end
names = {'SSA model vs data marginals', 'SSA model vs data: mean and std vs time'};
for i = 1:numel(names)
    close(findobj(groot, 'Type', 'figure', 'Name', names{i}));
end
M = model;
M.Solutions = state.sol;
M.ssaOptions.Nexp = 1;
M.ssaOptions.nSimsPerExpt = state.npe;
kdeArgs = ssit.parest.kdeArgsForRuns(kdeOpts, state.npe);
[~, figs] = M.computeLossFunctionSSA(lossFun, [], enforceIndependence, true, ...
    kdeArgs{:}, 'plotMarginals', true);
for i = 1:numel(figs)
    tl = findobj(figs(i), 'Type', 'tiledlayout');
    if ~isempty(tl)
        title(tl(1), sprintf('%s -- new best: log posterior estimate %.2f', ...
            tl(1).Title.String, value));
    end
end
if ~isempty(figurePrefix)
    suffix = {'_marginals', '_meanStd'};
    for i = 1:numel(figs)
        savefig(figs(i), [char(figurePrefix), suffix{i}, '.fig']);
        exportgraphics(figs(i), [char(figurePrefix), suffix{i}, '.png'], 'Resolution', 110);
    end
end
drawnow
end
