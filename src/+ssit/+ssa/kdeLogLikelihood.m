function [logL, pData, h, correction] = kdeLogLikelihood(samples, dataStates, dataCounts, opts)
% ssit.ssa.kdeLogLikelihood - Approximate log-likelihood of discrete
% single-cell data from SSA samples, using a kernel density estimate (KDE)
% to fill in states that the SSA samples did not visit.
%
%   [logL, pData, h] = ssit.ssa.kdeLogLikelihood(samples, dataStates, dataCounts)
%   [...] = ssit.ssa.kdeLogLikelihood(..., bandwidth=..., minBandwidth=..., probFloor=...)
%
% Inputs:
%   * samples    - (N x d) SSA samples (one row per trajectory, one column
%                  per observed species), nonnegative integers.
%   * dataStates - (M x d) observed states (nonnegative integers).
%   * dataCounts - (M x 1) number of cells observed in each state.
% Name-value options:
%   * bandwidth    - kernel widths (in molecule counts):
%       'lcv' (default)  likelihood cross-validation: h = c*hRef, with hRef
%                        the normal-reference widths below and the factor c
%                        chosen to maximize the leave-one-out log-likelihood
%                        of the SSA samples. Gives (almost) no smoothing
%                        when the samples are (nearly) deterministic and
%                        more smoothing where samples are sparse.
%       'silverman'      normal-reference rule
%                          h_j = s_j*(4/((d+2)*N))^(1/(d+4)),
%                        s_j = min(std, IQR/1.349) of species j.
%       numeric          scalar or 1 x d vector of widths.
%   * minBandwidth - lower limit on each h_j. Default: 0.5 for
%                    'silverman' (so constant species are still smoothed),
%                    0.01 for 'lcv' and numeric widths.
%   * probFloor    - (1e-10) lower limit on each probability before taking
%                    the log (the same floor used by SSIT.computeLikelihood).
%   * maxGrid      - (5e7) largest number of lattice points for the grid
%                    method; larger problems use the pairwise method.
%   * lcvFactors   - (logspace(-2,0.5,21)) candidate factors c for 'lcv'.
%   * biasCorrection - 'none' (default) or 'jensen'. With 'jensen', a
%                    delta-method correction for the downward bias of the
%                    log of a noisy estimate is added to each state's
%                    log-probability (see Method).
%
% Outputs:
%   * logL  - sum_m dataCounts(m)*(log(pData(m)) + correction(m)).
%   * pData - (M x 1) KDE probability of each observed state (after the
%             floor).
%   * h     - (1 x d) bandwidths used.
%   * correction - (M x 1) Jensen correction added to each log(pData(m))
%             (zeros when biasCorrection = 'none').
%
% Method:
%   The KDE is a mixture over the SSA samples x_k of product kernels
%     P(y) = (1/N) sum_k prod_j K_j(y_j | x_kj),
%   where K_j is a Gaussian of width h_j discretized onto the integers,
%     K_j(y | x) = Phi((y+0.5-x)/h_j) - Phi((y-0.5-x)/h_j)   for y >= 1,
%     K_j(0 | x) = Phi((0.5-x)/h_j)     (mass below 0 folded into 0),
%   so each kernel sums to one over the nonnegative integers and P is a
%   proper joint PMF. Because samples and data lie on the integer lattice,
%   P is computed as the sample histogram multiplied along each species by
%   a banded (sparse) kernel matrix, which is fast for the 1-3 observed
%   species typical of single-cell data. For very large lattices the
%   kernel sum is evaluated directly at the data states instead (pairwise
%   method; 'lcv' is then replaced by 'silverman').
%   The KDE is biased (smoothed); the bias shrinks as N grows and h -> 0.
%
%   Jensen correction ('jensen'): even where P is unbiased, log P is
%   biased low, E[log P(y)] ~ log E[P(y)] - Var(P(y))/(2 E[P(y)]^2), and
%   this term (always negative) is summed over every observed cell, which
%   makes it the largest part of the bias for small N. Because P(y) is a
%   mean of N i.i.d. kernel values K(y|x_k), its variance is estimated
%   without bias from the same samples,
%     Var(P(y)) ~ (P2(y) - P(y)^2)/(N-1),  P2(y) = (1/N) sum_k K(y|x_k)^2,
%   where P2 is computed like P but with squared kernels. The correction
%     correction(y) = Var(P(y))/(2 P(y)^2)
%   is added to log P(y). It lies in [0, 1/2] (since P2 <= N*P^2), so it can never
%   dominate the estimate. It is not applied to states at the probability
%   floor, nor when N < 2. It removes the noise part of the bias to leading
%   order but not the smoothing part.
arguments
    samples (:,:) double
    dataStates (:,:) double
    dataCounts (:,1) double
    opts.bandwidth = 'lcv'
    opts.minBandwidth = []
    opts.probFloor (1,1) double {mustBeNonnegative} = 1e-10
    opts.maxGrid (1,1) double {mustBePositive} = 5e7
    opts.lcvFactors (1,:) double {mustBePositive} = logspace(-2,0.5,21)
    opts.biasCorrection {mustBeMember(opts.biasCorrection,{'none','jensen'})} = 'none'
end

[N, d] = size(samples);
if N == 0
    error('SSIT:kdeLogLikelihood:noSamples', 'No SSA samples provided.');
end
if size(dataStates,2) ~= d
    error('SSIT:kdeLogLikelihood:dimMismatch', ...
        'samples has %d columns but dataStates has %d.', d, size(dataStates,2));
end
if size(dataStates,1) ~= numel(dataCounts)
    error('SSIT:kdeLogLikelihood:countMismatch', ...
        'dataStates has %d rows but dataCounts has %d entries.', ...
        size(dataStates,1), numel(dataCounts));
end
if any(samples(:) < 0) || any(dataStates(:) < 0)
    error('SSIT:kdeLogLikelihood:negative', 'States must be nonnegative integers.');
end

% Bandwidth rule and its default floor.
if ischar(opts.bandwidth) || isstring(opts.bandwidth)
    rule = lower(char(opts.bandwidth));
    if ~any(strcmp(rule, {'lcv','silverman'}))
        error('SSIT:kdeLogLikelihood:bandwidth', ...
            'Unknown bandwidth rule ''%s''. Use ''lcv'', ''silverman'', or numeric widths.', rule);
    end
else
    rule = 'numeric';
end
hMin = opts.minBandwidth;
if isempty(hMin)
    if strcmp(rule,'silverman'), hMin = 0.5; else, hMin = 0.01; end
end

% Normal-reference widths (used by 'silverman' and as the scale for 'lcv').
s = std(samples, 0, 1);
sortedSamples = sort(samples, 1);
q25 = sortedSamples(max(1, round(0.25*N)), :);
q75 = sortedSamples(max(1, round(0.75*N)), :);
iqrScale = (q75 - q25)/1.349;
useIqr = iqrScale > 0;
s(useIqr) = min(s(useIqr), iqrScale(useIqr));
hRef = s * (4/((d+2)*N))^(1/(d+4));

% Lattice for the grid method.
G = max([max(samples,[],1); max(dataStates,[],1)], [], 1) + 1;
useGrid = prod(G) <= opts.maxGrid;

switch rule
    case 'numeric'
        h = opts.bandwidth(:)';
        if isscalar(h)
            h = h*ones(1,d);
        elseif numel(h) ~= d
            error('SSIT:kdeLogLikelihood:bandwidth', ...
                'bandwidth must be a scalar or have one entry per species (%d).', d);
        end
        h = max(h, hMin);
    case 'silverman'
        h = max(hRef, hMin);
    case 'lcv'
        if useGrid
            h = lcvBandwidth(samples, G, hRef, hMin, opts.lcvFactors, opts.probFloor);
        else
            warning('SSIT:kdeLogLikelihood:lcvTooLarge', ...
                'Lattice too large for ''lcv''; using ''silverman'' widths.');
            h = max(hRef, max(hMin, 0.5));
        end
end
h = max(h, eps);

doJensen = strcmp(opts.biasCorrection, 'jensen') && N >= 2;
if useGrid
    if d == 1
        lin = dataStates + 1;
    else
        idx = num2cell(dataStates + 1, 1);
        lin = sub2ind(G, idx{:});
    end
    P = gridKDE(samples, G, h);
    pData = reshape(P(lin), [], 1);
    if doJensen
        P2 = gridKDE(samples, G, h, 2);
        p2Data = reshape(P2(lin), [], 1);
    end
else
    [pData, p2Data] = pairwiseKDE(samples, dataStates, h);
end

correction = zeros(size(pData));
if doJensen
    ok = pData > opts.probFloor;
    relVar = (p2Data(ok) - pData(ok).^2) ./ ((N-1) * pData(ok).^2);
    correction(ok) = min(max(relVar, 0), 1) / 2;
end

pData = max(pData, opts.probFloor);
logL = dataCounts' * (log(pData) + correction);
end

% -------------------------------------------------------------------------
function P = gridKDE(samples, G, h, power)
% KDE on the lattice 0..G-1: sample histogram multiplied along each
% dimension by the (sparse, banded) discretized-Gaussian kernel matrix.
% With power = 2, uses squared kernels, giving (1/N) sum_k K(y|x_k)^2.
if nargin < 4, power = 1; end
[N, d] = size(samples);
if d == 1
    P = accumarray(samples + 1, 1, [G 1])/N;
else
    P = accumarray(samples + 1, 1, G)/N;
end
for j = 1:d
    K = kernelMatrix(G(j), h(j));
    if power ~= 1, K = K.^power; end
    P = modeMultiply(P, K, j, G);
end
P = full(P);   % (sparse*full can stay sparse, e.g., for a 1-state lattice)
end

function K = kernelMatrix(G, h)
% Sparse G x G matrix K(y+1, x+1) = K(y | x) for y, x in 0..G-1, truncated
% where the kernel is negligible (|y-x| > 10*h + 1).
Phi = @(z) 0.5*erfc(-z/sqrt(2));
w = ceil(10*h) + 1;
x = 0:G-1;
[y, xx] = ndgrid(-w:w, x);
y = y + xx;                                   % target states near each source
keep = y >= 1 & y <= G-1;
y = y(keep); xx = xx(keep);
vals = Phi((y + 0.5 - xx)/h) - Phi((y - 0.5 - xx)/h);
% Row for y = 0 (all mass below 0.5).
x0 = x(x <= w + 1);
K = sparse([y; zeros(numel(x0),1)] + 1, [xx; x0(:)] + 1, ...
    [vals; Phi((0.5 - x0(:))/h)], G, G);
end

function P = modeMultiply(P, K, j, G)
% Multiply the j-th dimension of the d-dimensional array P by matrix K.
d = numel(G);
if d == 1
    P = K*P;
    return
end
order = [j, setdiff(1:d, j)];
Pj = reshape(permute(P, order), G(j), []);
Pj = K*Pj;
P = ipermute(reshape(Pj, G(order)), order);
end

function h = lcvBandwidth(samples, G, hRef, hMin, factors, probFloor)
% Choose h = c*hRef (floored at hMin) maximizing the leave-one-out
% log-likelihood of the samples:
%   sum_i log( (N*P(x_i) - Kself(x_i)) / (N-1) ),
% where Kself(x) = prod_j K_j(x_j | x_j) is the sample's own contribution.
[N, d] = size(samples);
if N < 2
    h = max(hRef, hMin);
    return
end
[uS, ~, ic] = unique(samples, 'rows');
counts = accumarray(ic, 1);
if d == 1
    lin = uS + 1;
else
    idx = num2cell(uS + 1, 1);
    lin = sub2ind(G, idx{:});
end
Phi = @(z) 0.5*erfc(-z/sqrt(2));
best = -Inf;
h = max(hRef, hMin);
for c = factors
    hc = max(c*hRef, hMin);
    P = gridKDE(samples, G, hc);
    kSelf = ones(size(uS,1), 1);
    for j = 1:d
        kj = Phi(0.5/hc(j)) - Phi(-0.5/hc(j))*ones(size(uS,1),1);
        kj(uS(:,j) == 0) = Phi(0.5/hc(j));   % boundary state
        kSelf = kSelf .* kj;
    end
    pLoo = (N*P(lin) - kSelf)/(N - 1);
    score = counts' * log(max(pLoo, probFloor));
    if score > best
        best = score;
        h = hc;
    end
end
end

function [pData, p2Data] = pairwiseKDE(samples, dataStates, h)
% Direct evaluation at the data states (for lattices too large to grid).
[N, d] = size(samples);
[uS, ~, ic] = unique(samples, 'rows');
w = accumarray(ic, 1)/N;
K = size(uS, 1);
M = size(dataStates, 1);
pData = zeros(M, 1);
p2Data = zeros(M, 1);
chunk = max(1, floor(2e7/K));
Phi = @(z) 0.5*erfc(-z/sqrt(2));
for i0 = 1:chunk:M
    rows = i0:min(M, i0+chunk-1);
    kern = ones(numel(rows), K);
    for j = 1:d
        y = dataStates(rows, j);
        x = uS(:, j)';
        upper = Phi((y + 0.5 - x)/h(j));
        lower = Phi((y - 0.5 - x)/h(j));
        lower(y <= 0, :) = 0;
        kern = kern .* (upper - lower);
    end
    pData(rows) = kern*w;
    p2Data(rows) = (kern.^2)*w;
end
end
