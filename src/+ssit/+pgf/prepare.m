function prob = prepare(src, opts, tSpan, model)
% prob = ssit.pgf.prepare(src, opts, tSpan) translates an SSIT model for the
% PGF solver and resolves the options that do not depend on the grid.
% prob = ssit.pgf.prepare(src, opts, tSpan, model) reuses a model already
% translated by ssit.pgf.translateModel (same structure; only the
% parameter values are taken from src).
% Shared by ssit.pgf.solve and the fitting routines (ssit.pgf.fitSetup),
% which then call ssit.pgf.chooseGrid / ssit.pgf.setGrid and
% ssit.pgf.propagate.
%
%   src   : SSIT object or struct (see ssit.pgf.solve).
%   opts  : PGF options ([] -> src.pgfOptions / ssit.pgf.defaultOptions).
%   tSpan : output times ([] -> src.tSpan).
%
% Returns struct prob with fields
%   model       : ssit.pgf.translateModel result,
%   theta       : parameter values (column, order of src.parameters),
%   mu          : 1 x dY Smolyak layers of the PGF species,
%   lowerFixed  : 1 x dY lower domain ends, or [] for lowerBound = 'auto',
%   init        : struct x0, p0, steadyState,
%   opts, tSpan,
%   sg, ops, lower : the grid, operators and lower ends (empty until
%                 ssit.pgf.setGrid is called).
arguments
    src
    opts = []
    tSpan = []
    model = []
end
opts = ssit.pgf.mergeOptions(src, opts);
if isempty(tSpan)
    tSpan = src.tSpan;
end
tSpan = tSpan(:)';
if isempty(model)
    model = ssit.pgf.translateModel(src, tSpan, opts);
else
    model.parameterValues = [src.parameters{:,2}]';
end
d = model.d;
Y = model.pgfIdx;

mu = perPgfSpecies(opts.mu, d, Y, 'mu');
if ischar(opts.lowerBound) || isstring(opts.lowerBound)
    if ~strcmpi(opts.lowerBound, 'auto')
        error('ssit:pgf:options', 'pgfOptions.lowerBound must be ''auto'' or numeric.');
    end
    lowerFixed = [];
else
    lowerFixed = perPgfSpecies(opts.lowerBound, d, Y, 'lowerBound');
end

x0 = src.initialCondition;
p0 = 1;
if isPropOrField(src, 'initialProbs') && ~isempty(src.initialProbs)
    p0 = src.initialProbs;
end
useSS = opts.initialSteadyState;
if isempty(useSS)
    useSS = isPropOrField(src, 'fspOptions') && isfield(src.fspOptions, 'initApproxSS') && ...
        ~isempty(src.fspOptions.initApproxSS) && src.fspOptions.initApproxSS;
end
init = struct('x0', x0, 'p0', p0, 'steadyState', logical(useSS));

prob = struct('model', model, 'theta', model.parameterValues(:), 'mu', mu, ...
    'lowerFixed', lowerFixed, 'init', init, 'opts', opts, 'tSpan', tSpan, ...
    'sg', [], 'ops', [], 'lower', []);
end

function v = perPgfSpecies(v, d, Y, name)
% Expand a scalar / per-species / per-PGF-species option to one entry per PGF species.
v = v(:)';
dY = numel(Y);
if isscalar(v)
    v = repmat(v, 1, dY);
elseif numel(v) == d && d ~= dY
    v = v(Y);
elseif numel(v) ~= dY
    error('ssit:pgf:options', ...
        'pgfOptions.%s must be a scalar or have one entry per PGF species (%d) or per species (%d).', ...
        name, dY, d);
end
end

function tf = isPropOrField(obj, name)
tf = (isstruct(obj) && isfield(obj, name)) || (isobject(obj) && isprop(obj, name));
end
