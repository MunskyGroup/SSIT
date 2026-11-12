function fp = fitSetup(src, opts)
% fp = ssit.pgf.fitSetup(src) prepares the PGF fit of an SSIT model to the
% data loaded with SSIT.loadData.  Everything that does not depend on the
% parameters is done once here: model translation, grid and operators
% (domain fixed for the whole fit), empirical PGF of the data and the
% residual weights.  Evaluate with ssit.pgf.residuals / ssit.pgf.logLikelihood,
% optimise with ssit.pgf.levenbergMarquardt or SSIT.maximizeLikelihood.
%
%   src  : SSIT object with loaded data (dataSet), fittingOptions
%          (modelVarsToFit, timesToFit, logPrior), initialTime.
%   opts : PGF options ([] -> src.pgfOptions); fields loss, errorFloor,
%          mu, lowerBound, ... (see ssit.pgf.defaultOptions).
%
% What is compared.  For a cell with state x, the data provide the
% empirical PGF of the observed species
%   Hd_g(z) = (1/n) sum_cells 1{observed enumerated values = g} prod_j z_j^x_j
% (one group g per combination of values of the observed enumerated
% species; a single group if none is observed).  The model counterpart is
% the sum of the blocks G_e over the enumerated states in group g, with the
% unobserved PGF species at z = 1 (through the interpolant).  The points z
% are the grid nodes (with unobserved coordinates set to 1).  Unobserved
% enumerated species are summed over.
%
% Loss (opts.loss), r = model - data at one time:
%   'gls'      logL = -1/2 sum_t r' (C_t/n_t + s^2 I)^-1 r, where C_t is the
%              covariance of the vector (1{g} z^x) over cells (from the
%              data) and s = errorFloor.  Gaussian (CLT) approximation of the
%              likelihood of the empirical PGF: the default.
%   'weighted' as 'gls' with the diagonal of C_t only.
%   'mse'      logL = -mean of r.^2 over times and points (Edward's objective).
% A log prior (fittingOptions.logPrior) is added to logL.
%
% Returns struct fp with fields
%   prob        : solver problem (ssit.pgf.prepare, grid set),
%   theta0      : all parameter values,  fitIdx : indices of fitted ones,
%   parameterNames, observedSpecies, fitTimes, tCol (columns of prob.tSpan),
%   Phi, Agg, M, nG, groupValues : comparison points and enumerated groups,
%   data(k).m, data(k).W, nCells : empirical PGF and weights per time,
%   loss, scale (logL = -scale*sum(r.^2)), errorFloor, logPrior, nRes.
arguments
    src
    opts = []
end
opts = ssit.pgf.mergeOptions(src, opts);
loss = lower(char(opts.loss));
if ~ismember(loss, {'gls', 'weighted', 'mse'})
    error('ssit:pgf:options', 'pgfOptions.loss must be ''gls'', ''weighted'' or ''mse''.');
end

% ---- data ---------------------------------------------------------------
species = src.species(:)';
ds = src.dataSet;
if isempty(ds) || ~isfield(ds, 'app') || ~isfield(ds, 'times')
    error('ssit:pgf:noData', 'No data loaded: call SSIT.loadData before fitting.');
end
T = ds.app.DataLoadingAndFittingTabOutputs.dataTensor;
subs = double(T.subs);
vals = double(T.vals);
times = ds.times(:)';
linked = linkedSpeciesNames(ds);
if isfield(ds, 'tensorSpecies') && ~isempty(ds.tensorSpecies)
    tensorSpecies = cellstr(ds.tensorSpecies);
else
    tensorSpecies = species(ismember(species, linked));
end
if size(subs, 2) - 1 ~= numel(tensorSpecies)
    error('ssit:pgf:data', ...
        'Cannot match the %d data columns to model species; reload the data with SSIT.loadData.', ...
        size(subs, 2) - 1);
end
obsDims = find(ismember(tensorSpecies, linked) & ismember(tensorSpecies, species));
if isempty(obsDims)
    error('ssit:pgf:data', 'None of the model species is linked to a data column.');
end
[~, obsSpecies] = ismember(tensorSpecies(obsDims), species);

if isprop(src, 'pdoOptions') || (isstruct(src) && isfield(src, 'pdoOptions'))
    if isfield(src.pdoOptions, 'PDO') && ~isempty(src.pdoOptions.PDO)
        error('ssit:pgf:pdo', 'Probabilistic distortion operators (pdoOptions.PDO) are not supported by the PGF fit.');
    end
end

% ---- times and parameters -------------------------------------------------
fo = src.fittingOptions;
tIdx = indexList(fo.timesToFit, numel(times));
fitTimes = times(tIdx);
t0 = src.initialTime;
if any(fitTimes < t0 - 1e-12 * max(1, abs(t0)))
    error('ssit:pgf:data', ['Data times (%s) precede initialTime = %g. Exclude them with ', ...
        'fittingOptions.timesToFit or change initialTime.'], num2str(fitTimes(fitTimes < t0)), t0);
end
tSpan = unique([t0, fitTimes]);
[~, tCol] = ismember(fitTimes, tSpan);
fitIdx = indexList(fo.modelVarsToFit, size(src.parameters, 1));

prob = ssit.pgf.prepare(src, opts, tSpan);
model = prob.model;
theta0 = prob.theta;
Y = model.pgfIdx;
E = model.enumIdx;
dY = numel(Y);
nE = size(model.enumStates, 2);
[obsY, dimY] = ismember(Y, obsSpecies);       % dimY: column of the data counts
[obsE, dimE] = ismember(E, obsSpecies);

nF = numel(tIdx);
X = cell(1, nF);
w = cell(1, nF);
nCells = zeros(1, nF);
for k = 1:nF
    rows = subs(:,1) == tIdx(k);
    X{k} = subs(rows, obsDims + 1) - 1;        % counts, columns = obsSpecies
    w{k} = vals(rows);
    nCells(k) = sum(w{k});
    if nCells(k) <= 0
        error('ssit:pgf:data', 'No data at time %g.', fitTimes(k));
    end
end

% ---- domain (fixed during the fit) --------------------------------------
if isempty(prob.lowerFixed)
    lowerB = 0.9 * ones(1, dY);
    if any(~obsY)
        probA = ssit.pgf.chooseGrid(prob, theta0);      % unobserved species: model means
        lowerB = probA.lower;
    end
    for j = find(obsY)
        m = 0;
        for k = 1:nF
            m = max(m, (w{k}.' * X{k}(:, dimY(j))) / nCells(k));
        end
        lowerB(j) = min(max(-1, 1 - opts.autoScale / max(m, eps)), 0.999);
    end
    lowerB(model.secondOrder(Y)) = -1;
    prob = ssit.pgf.setGrid(prob, lowerB);
else
    prob = ssit.pgf.setGrid(prob, prob.lowerFixed);
end
G0 = ssit.pgf.propagate(prob, theta0);          % also checks stability once
sg = prob.sg;
N = sg.N;

% ---- comparison points and enumerated groups -------------------------------
Z = sg.z;
Z(:, ~obsY) = 1;
if all(obsY)
    Phi = [];
    Zu = Z;
else
    Zu = unique(Z, 'rows', 'stable');
    Phi = ssit.pgf.smolyakBasis(Zu, sg.multiIndex, sg.domain) * sg.invInterp;
end
M = size(Zu, 1);
Zobs = Zu(:, obsY);
colY = dimY(obsY);
if any(obsE)
    colE = dimE(obsE);
    modelKeys = model.enumStates(obsE, :).';
    dataKeys = cell2mat(cellfun(@(x) x(:, colE), X(:), 'UniformOutput', false));
    groupValues = unique([modelKeys; dataKeys], 'rows');
    [~, gModel] = ismember(modelKeys, groupValues, 'rows');
    nG = size(groupValues, 1);
    Agg = full(sparse(gModel(:), (1:nE)', 1, nG, nE));
else
    colE = [];
    groupValues = zeros(1, 0);
    nG = 1;
    Agg = ones(1, nE);
end

% ---- empirical PGF and covariances ---------------------------------------
data = struct('m', cell(1, nF), 'C', [], 'W', []);
for k = 1:nF
    B = ones(size(X{k}, 1), M);
    for jj = 1:numel(colY)
        B = B .* ((Zobs(:,jj).') .^ X{k}(:, colY(jj)));
    end
    if nG > 1
        [~, g] = ismember(X{k}(:, colE), groupValues, 'rows');
        Yr = zeros(size(B, 1), M*nG);
        for gg = 1:nG
            Yr(g == gg, (gg-1)*M + (1:M)) = B(g == gg, :);
        end
    else
        Yr = B;
    end
    n = nCells(k);
    m = (w{k}.' * Yr).' / n;
    data(k).m = m;
    switch loss
        case 'gls'
            data(k).C = (Yr.' * (Yr .* w{k})) / n - m * m.';
        case 'weighted'
            data(k).C = (w{k}.' * (Yr.^2)).' / n - m.^2;
    end
end

% ---- error floor and weights ------------------------------------------------
s = 0;
if ~strcmp(loss, 'mse')
    s = opts.errorFloor;
    if ischar(s) || isstring(s)
        s = estimateErrorFloor(prob, theta0, G0, tCol, Agg, Phi, Zu, nE, opts);
    elseif ~(isscalar(s) && s > 0)
        error('ssit:pgf:options', 'pgfOptions.errorFloor must be ''auto'' or a positive number.');
    end
end
for k = 1:nF
    n = nCells(k);
    switch loss
        case 'gls'
            C = (data(k).C + data(k).C.') / (2*n);
            [V, D] = eig(C);
            lam = max(diag(D), 0);
            data(k).W = (V ./ sqrt(lam + s^2).').';     % W'W = (C + s^2 I)^-1
        case 'weighted'
            data(k).W = 1 ./ sqrt(max(data(k).C, 0) / n + s^2);
        case 'mse'
            data(k).W = 1 / sqrt(nF * M * nG);
    end
    data(k).C = [];
end

logPrior = [];
if isfield(fo, 'logPrior')
    logPrior = fo.logPrior;
end
nPrior = 0;
if ~isempty(logPrior)
    nPrior = numel(logPrior(theta0(fitIdx)));
end
scale = 0.5;
if strcmp(loss, 'mse')
    scale = 1;
end
fp = struct('prob', prob, 'theta0', theta0, 'fitIdx', fitIdx(:)', ...
    'parameterNames', {src.parameters(fitIdx, 1)'}, 'observedSpecies', {species(obsSpecies)}, ...
    'fitTimes', fitTimes, 'tCol', tCol, 'Phi', Phi, 'Agg', Agg, 'M', M, 'nG', nG, ...
    'groupValues', groupValues, 'points', Zu, 'data', {data}, 'nCells', nCells, 'loss', loss, ...
    'scale', scale, 'errorFloor', s, 'logPrior', logPrior, 'nRes', nF*M*nG + nPrior);
if opts.verbose
    fprintf(['PGF fit setup: %d times, observed {%s}, enumerated {%s} (%d groups), N = %d nodes, ', ...
        '%d points, lowerBound = [%s], loss = %s, errorFloor = %.3g\n'], nF, ...
        strjoin(fp.observedSpecies, ', '), strjoin(model.species(E), ', '), nG, N, M, ...
        num2str(prob.lower, '%.4g '), loss, s);
end
end

% -------------------------------------------------------------------------
function s = estimateErrorFloor(prob, theta0, G0, tCol, Agg, Phi, Zu, nE, opts)
% RMS change of the compared PGF values when every mu is raised by one.
try
    probF = prob;
    probF.mu = prob.mu + 1;
    probF = ssit.pgf.setGrid(probF, prob.lower);
    GF = ssit.pgf.propagate(probF, theta0, false);
    PhiF = ssit.pgf.smolyakBasis(Zu, probF.sg.multiIndex, probF.sg.domain) * probF.sg.invInterp;
    Hc = ssit.pgf.comparisonValues(G0, prob.sg.N, nE, tCol, Agg, Phi);
    Hf = ssit.pgf.comparisonValues(GF, probF.sg.N, nE, tCol, Agg, PhiF);
    s = sqrt(mean((Hf(:) - Hc(:)).^2));
catch err
    s = 1e-4;
    warning('ssit:pgf:errorFloor', ...
        'Could not estimate the PGF error floor on a refined grid (%s); using %g.', err.message, s);
end
s = max(s, opts.relTol);
end

function idx = indexList(v, n)
% 'all' / logical mask / index list -> row of indices.
if ischar(v) || isstring(v)
    if strcmpi(v, 'all')
        idx = 1:n;
        return
    end
    error('ssit:pgf:options', 'Expected ''all'', a logical mask or indices.');
end
if islogical(v)
    idx = find(v(:)');
else
    idx = v(:)';
end
end

function names = linkedSpeciesNames(ds)
names = {};
if ~isfield(ds, 'linkedSpecies') || isempty(ds.linkedSpecies)
    return
end
ls = ds.linkedSpecies;
if isa(ls, 'dictionary')
    names = cellstr(keys(ls))';
elseif iscell(ls)
    names = ls(:,1)';
end
names = cellfun(@char, names, 'UniformOutput', false);
end
