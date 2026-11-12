function [prob, G] = chooseGrid(prob, theta)
% [prob, G] = ssit.pgf.chooseGrid(prob, theta) sets the z-domain of the PGF
% species and returns the block PGF G (nE*N x Nt) solved on it at theta
% (default prob.theta).
%   lowerBound numeric : that domain (validated by ssit.pgf.setGrid).
%   lowerBound 'auto'  : solve on a provisional domain, then
%                        lower_j = max(-1, 1 - autoScale / max_t E[y_j]),
%                        capped at 0.999, and re-solve (autoPasses times).
%                        Species with quadratic propensity terms get -1.
% If prob already holds a grid (prob.sg, e.g. from a previous solve, see
% ssit.pgf.solve), it is used when it matches the requested domain:
%   numeric : reused if prob.lower equals lowerBound,
%   'auto'  : used as the provisional domain, and kept without re-solving
%             if the domain it implies is close to it (1 - lower within a
%             factor 1.5 for every species).
arguments
    prob (1,1) struct
    theta = []
end
if isempty(theta)
    theta = prob.theta;
end
haveGrid = ~isempty(prob.sg);
if ~isempty(prob.lowerFixed)
    if ~(haveGrid && isequal(prob.lower, prob.lowerFixed))
        prob = ssit.pgf.setGrid(prob, prob.lowerFixed);
    end
    G = ssit.pgf.propagate(prob, theta);
    return
end
opts = prob.opts;
model = prob.model;
Y = model.pgfIdx;
secondY = model.secondOrder(Y);
if ~haveGrid
    lower = 0.9 * ones(1, numel(Y));
    lower(secondY) = -1;
    prob = ssit.pgf.setGrid(prob, lower);
end
G = ssit.pgf.propagate(prob, theta);
for pass = 1:opts.autoPasses
    means = ssit.pgf.computeMoments(model, prob.sg, G);
    lower = autoLower(max(means(Y,:), [], 2)', opts.autoScale);
    lower(secondY) = -1;
    if haveGrid && pass == 1 && domainFits(lower, prob.lower)
        return
    end
    prob = ssit.pgf.setGrid(prob, lower);
    G = ssit.pgf.propagate(prob, theta);
end
end

function lower = autoLower(m, autoScale)
m = max(m, 0);
lower = max(-1, 1 - autoScale ./ max(m, eps));
lower = min(lower, 0.999);
end

function tf = domainFits(target, current)
% Widths 1 - lower agree within a factor 1.5 for every species.
ratio = (1 - target) ./ (1 - current);
tf = all(ratio >= 1/1.5 & ratio <= 1.5);
end
