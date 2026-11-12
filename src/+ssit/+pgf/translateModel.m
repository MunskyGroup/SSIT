function model = translateModel(src, tSpan, opts)
% model = ssit.pgf.translateModel(src, tSpan, opts) converts an SSIT reaction
% model into the reaction terms of the (block) PGF equation.
%
% Species are split into ENUMERATED species E (small, bounded counts such as
% gene states; their reachable states e are listed explicitly) and PGF
% species Y.  For each enumerated state e the PGF of Y restricted to e is
%   G_e(z) = E[ 1{x_E = e} prod_j z_j^y_j ],
% and a reaction r with stoichiometry (s_E, s_Y), whose propensity at fixed e
% is written in falling factorials of y,
%   w_r(e, y) = sum_k c_{r,k}(e, t, theta) prod_j y_j^(k_j),
% contributes
%   dG_{e+s_E}/dt += c_{r,k}(e) z^(s_Y + k) d^k G_e,
%   dG_e/dt       -= c_{r,k}(e) z^k d^k G_e.
% With no enumerated species there is a single block (phase-2 behaviour).
%
%   src   : SSIT object (or struct): species, propensityFunctions,
%           stoichiometry, parameters, inputExpressions, initialCondition.
%   tSpan : times used to probe time-dependent propensities.
%   opts  : PGF options (enumeratedSpecies, maxEnumeratedCount,
%           maxEnumeratedStates); default src.pgfOptions / defaultOptions.
%
% The translation is numeric (no Symbolic Toolbox): for a polynomial of
% degree <= K_j in each PGF species, c_k = Delta^k w(0)/k! exactly, a fixed
% linear combination of propensity values on the lattice {0..K_1} x ... .
% Propensities may be any function of time and of the enumerated species.
% The structure is probed with randomised parameter values, so it does not
% depend on the current values (e.g. a rate set to 0).
%
% Errors
%   ssit:pgf:nonPolynomial      propensity not polynomial (degree <= 2) in a PGF species
%   ssit:pgf:negativeState      propensity nonzero where a count would go negative
%   ssit:pgf:enumeratedUnbounded  enumerated species reach more than maxEnumeratedStates states
%   ssit:pgf:noPgfSpecies       every species is enumerated (use FSP)
%
% Returns struct model with fields
%   species, d, R, S, parameterNames, parameterValues,
%   enumIdx, pgfIdx        : indices of enumerated / PGF species,
%   enumStates             : dE x nE reachable enumerated states (nE = 1, dE = 0 if none),
%   rxn(r).w               : @(x,t,k) vectorised propensity (x is d x n),
%   rxn(r).X               : dY x nj lattice over the PGF species,
%   rxn(r).Xfull           : d x (nj*nE) evaluation points (lattice x enumerated states),
%   rxn(r).k, rxn(r).A     : falling-factorial orders kept (nT x dY) and weights (nT x nj),
%   rxn(r).target          : 1 x nE index of e + s_E (0 where the reaction cannot fire),
%   rxn(r).timeDependent,
%   secondOrder, mixedOrder: 1 x d flags (PGF species in quadratic / mixed terms),
%   timeDependent.
arguments
    src
    tSpan (1,:) double = [0 1]
    opts = []
end
opts = mergeOptions(src, opts);
species = src.species(:)';
d = numel(species);
S = src.stoichiometry;
if size(S,1) ~= d
    error('ssit:pgf:stoichiometry', 'The stoichiometry must have one row per species.');
end
R = size(S,2);
props = src.propensityFunctions(:);
if numel(props) ~= R
    error('ssit:pgf:stoichiometry', 'There are %d propensities but %d stoichiometry columns.', numel(props), R);
end
parNames = src.parameters(:,1)';
parValues = [src.parameters{:,2}]';
inputs = {};
if isPropOrField(src, 'inputExpressions')
    inputs = src.inputExpressions;
end
x0 = src.initialCondition;
if size(x0,1) ~= d
    error('ssit:pgf:initialCondition', 'initialCondition must have one row per species.');
end

w = cell(1,R);
for r = 1:R
    w{r} = makePropensityFunction(props{r}, species, parNames, inputs);
end

% Probe settings (private random stream: the global RNG is not touched).
pr.rs = RandStream('mt19937ar', 'Seed', 20261106);
pr.pars = {probeParameters(parValues, pr.rs), probeParameters(parValues, pr.rs)};
t0 = tSpan(1); t1 = tSpan(end);
pr.times = unique([tSpan, linspace(t0, t1, 7), t0 + (t1 - t0) * rand(pr.rs, 1, 4)]);
pr.nBase = 6;  pr.maxBase = 25;  pr.nValidate = 60;  pr.maxValidate = 60;
pr.tol = 1e-8;

% ---- enumerated species and their reachable states ----
E = chooseEnumerated(opts, species, S, x0, w, props, pr);
enumIdx = find(E);
pgfIdx = find(~E);
if isempty(pgfIdx)
    error('ssit:pgf:noPgfSpecies', ...
        ['All species are enumerated (%s): the model is a finite Markov chain, ', ...
         'use the FSP solver.'], strjoin(species, ', '));
end
[enumStates, fires] = enumerateStates(enumIdx, d, S, x0, w, props, species, opts.maxEnumeratedStates, pr);
nE = size(enumStates, 2);
dY = numel(pgfIdx);

model = struct('species', {species}, 'd', d, 'R', R, 'S', S, ...
    'parameterNames', {parNames}, 'parameterValues', parValues, ...
    'enumIdx', enumIdx, 'pgfIdx', pgfIdx, 'enumStates', enumStates, ...
    'rxn', [], 'secondOrder', false(1,d), 'mixedOrder', false(1,d), 'timeDependent', false);
rxn = repmat(struct('w', [], 'expression', '', 'X', [], 'Xfull', [], 'k', [], 'A', [], ...
    'target', [], 'timeDependent', false), 1, R);

for r = 1:R
    expr = props{r};
    rxn(r).w = w{r};
    rxn(r).expression = expr;

    % --- degree in each PGF species, for every enumerated state ---
    K = zeros(1,dY);
    for jj = 1:dY
        i = pgfIdx(jj);
        for ip = 1:numel(pr.pars)
            for t = pr.times
                B = randomStates(d, pr.nBase, pr.maxBase, enumIdx, enumStates, pr.rs);
                n = size(B,2);
                V = zeros(4, n);
                for step = 0:3
                    Bs = B; Bs(i,:) = Bs(i,:) + step;
                    V(step+1,:) = evalW(w{r}, Bs, t, pr.pars{ip}, r, expr);
                end
                scale = max(abs(V), [], 'all') + realmin;
                D1 = V(2,:) - V(1,:);
                D2 = V(3,:) - 2*V(2,:) + V(1,:);
                D3 = V(4,:) - 3*V(3,:) + 3*V(2,:) - V(1,:);
                if any(abs(D3) > pr.tol*scale)
                    error('ssit:pgf:nonPolynomial', ...
                        ['Reaction %d (propensity ''%s'') is not a polynomial of degree <= 2 in species ''%s''. ', ...
                         'The PGF solver needs propensities that are polynomial in the non-enumerated ', ...
                         'species (e.g. mass action); enumerate ''%s'' if its count is small ', ...
                         '(pgfOptions.enumeratedSpecies), or use the FSP or SSA solvers.'], ...
                        r, expr, species{i}, species{i});
                end
                if any(abs(D2) > pr.tol*scale)
                    K(jj) = max(K(jj), 2);
                elseif any(abs(D1) > pr.tol*scale)
                    K(jj) = max(K(jj), 1);
                end
            end
        end
    end

    % --- lattice over the PGF species and falling-factorial weights ---
    X = latticePoints(K);                       % dY x nj
    nj = size(X,2);
    kAll = X';
    A = zeros(size(kAll,1), nj);
    for q = 1:size(kAll,1)
        kq = kAll(q,:);
        for j = 1:nj
            xj = X(:,j)';
            if all(xj <= kq)
                A(q,j) = prod((-1).^(kq - xj) .* arrayfun(@nchoosek, kq, xj) ./ factorial(kq));
            end
        end
    end
    Xfull = zeros(d, nj*nE);
    Xfull(pgfIdx,:) = repmat(X, 1, nE);
    Xfull(enumIdx,:) = kron(enumStates, ones(1,nj));

    % --- coefficients at the probes (terms x states x probes x times) ---
    nP = numel(pr.pars); nT = numel(pr.times);
    C = zeros(size(kAll,1), nE, nP, nT);
    for ip = 1:nP
        for it = 1:nT
            v = evalW(w{r}, Xfull, pr.times(it), pr.pars{ip}, r, expr);
            C(:,:,ip,it) = A * reshape(v, nj, nE);
        end
    end
    cScale = max(abs(C), [], 'all') + realmin;
    nonzero = max(abs(C), [], [3 4]) > pr.tol*cScale;     % terms x states
    keep = any(nonzero, 2);                               % identically zero -> no terms
    for ip = 1:nP
        Cp = reshape(C(keep,:,ip,:), nnz(keep)*nE, nT);
        if any(abs(Cp - Cp(:,1)) > pr.tol*cScale, 'all')
            rxn(r).timeDependent = true;
        end
    end
    kKept = kAll(keep,:);
    A = A(keep,:);
    nonzero = nonzero(keep,:);

    % --- validation at random states: w must equal its expansion in every state ---
    for ip = 1:nP
        for t = pr.times
            Xv = randomStates(d, pr.nValidate, pr.maxValidate, enumIdx, enumStates, pr.rs);
            wv = evalW(w{r}, Xv, t, pr.pars{ip}, r, expr);
            c = A * reshape(evalW(w{r}, Xfull, t, pr.pars{ip}, r, expr), nj, nE);   % nT x nE
            if isempty(enumIdx)
                eOf = ones(1, size(Xv,2));
            else
                [~, eOf] = ismember(Xv(enumIdx,:)', enumStates', 'rows');
                eOf = eOf';
            end
            wp = sum(c(:, eOf) .* fallingFactorials(Xv(pgfIdx,:), kKept), 1);
            if any(abs(wv - wp) > pr.tol*(max(abs(wv)) + realmin))
                error('ssit:pgf:nonPolynomial', ...
                    ['Reaction %d (propensity ''%s'') is not a polynomial in the non-enumerated species. ', ...
                     'The PGF solver needs polynomial propensities there; use FSP or SSA.'], r, expr);
            end
        end
    end

    % --- no negative counts of PGF species ---
    sY = S(pgfIdx, r)';
    for q = 1:size(kKept,1)
        bad = find(sY + kKept(q,:) < 0, 1);
        if ~isempty(bad)
            sp = species{pgfIdx(bad)};
            error('ssit:pgf:negativeState', ...
                ['Reaction %d (propensity ''%s'') is nonzero when species ''%s'' has fewer than %d molecules, ', ...
                 'so it would create negative counts. Write the propensity so that it vanishes there ', ...
                 '(e.g. ''g*%s'' instead of ''g'').'], r, expr, sp, -sY(bad), sp);
        end
    end

    % --- target enumerated state of each firing ---
    target = zeros(1, nE);
    sE = S(enumIdx, r);
    for e = 1:nE
        if fires(r,e) || any(nonzero(:,e))
            if isempty(enumIdx)
                target(e) = 1;
            else
                [~, target(e)] = ismember((enumStates(:,e) + sE)', enumStates', 'rows');
                if target(e) == 0
                    error('ssit:pgf:negativeState', ...
                        'Reaction %d (propensity ''%s'') leads outside the reachable enumerated states.', r, expr);
                end
            end
        end
    end

    rxn(r).X = X;
    rxn(r).Xfull = Xfull;
    rxn(r).k = kKept;
    rxn(r).A = A;
    rxn(r).target = target;
    sec = any(kKept >= 2, 1);
    mix = any(kKept(sum(kKept,2) >= 2 & max(kKept,[],2) < 2, :) > 0, 1);
    model.secondOrder(pgfIdx) = model.secondOrder(pgfIdx) | sec;
    model.mixedOrder(pgfIdx) = model.mixedOrder(pgfIdx) | mix;
end
model.rxn = rxn;
model.timeDependent = any([rxn.timeDependent]);
end

% =========================================================================
function E = chooseEnumerated(opts, species, S, x0, w, props, pr)
d = numel(species);
spec = opts.enumeratedSpecies;
if ischar(spec) || isstring(spec)
    if ~strcmpi(spec, 'auto')
        spec = cellstr(spec);
    end
end
if iscell(spec) || (isstring(spec) && ~strcmpi(spec, 'auto'))
    spec = cellstr(spec);
    [found, loc] = ismember(spec, species);
    if ~all(found)
        error('ssit:pgf:options', 'Unknown species in pgfOptions.enumeratedSpecies: %s.', ...
            strjoin(spec(~found), ', '));
    end
    E = false(1,d);
    E(loc) = true;
    return
end

% ---- automatic detection ----
maxCount = opts.maxEnumeratedCount;
E = false(1,d);
% (a) nonnegative conservation laws c'S = 0 with a small total c'x0.
Nul = null(S', 'r');
cands = Nul;
m = size(Nul,2);
for a = 1:m
    for b = a+1:m
        cands = [cands, Nul(:,a) + Nul(:,b), Nul(:,a) - Nul(:,b)]; %#ok<AGROW>
    end
end
cands = [cands, -cands];
for c = cands
    if all(c >= -1e-12) && any(c > 1e-12)
        c(c < 0) = 0;
        T = max(c' * x0);
        bound = floor(T ./ c(c > 1e-12) + 1e-9);
        idx = find(c > 1e-12);
        E(idx(bound <= maxCount)) = true;
    end
end
% (b) species bounded by their own propensities: every reaction that
%     increases species i vanishes once x_i is within s_ir of a bound U.
for i = find(~E)
    xmax = max(x0(i,:));
    if xmax > maxCount
        continue
    end
    inc = find(S(i,:) > 0);
    if isempty(inc)
        E(i) = true;          % never increases: bounded by its initial count
        continue
    end
    for U = xmax:maxCount
        ok = true;
        for r = inc
            for xi = max(0, U - S(i,r) + 1):U
                for ip = 1:numel(pr.pars)
                    for t = pr.times
                        B = randi(pr.rs, [0 pr.maxBase], d, pr.nBase);
                        B(i,:) = xi;
                        if any(evalW(w{r}, B, t, pr.pars{ip}, r, props{r}) ~= 0)
                            ok = false; break
                        end
                    end
                    if ~ok, break, end
                end
                if ~ok, break, end
            end
            if ~ok, break, end
        end
        if ok
            E(i) = true;
            break
        end
    end
end
end

function [states, fires] = enumerateStates(enumIdx, d, S, x0, w, props, species, maxStates, pr)
% Breadth-first search of the enumerated states reachable from the initial
% condition.  fires(r,e) = reaction r can have a nonzero propensity in state e.
R = size(S,2);
if isempty(enumIdx)
    states = zeros(0,1);
    fires = true(R,1);
    return
end
states = unique(x0(enumIdx,:)', 'rows', 'stable')';
fires = false(R, 0);
head = 1;
while head <= size(states,2)
    e = states(:,head);
    for r = 1:R
        can = false;
        for ip = 1:numel(pr.pars)
            for t = pr.times
                B = randi(pr.rs, [0 pr.maxBase], d, pr.nBase);
                B(enumIdx,:) = repmat(e, 1, pr.nBase);
                if any(evalW(w{r}, B, t, pr.pars{ip}, r, props{r}) ~= 0)
                    can = true; break
                end
            end
            if can, break, end
        end
        fires(r,head) = can;
        if ~can
            continue
        end
        tgt = e + S(enumIdx, r);
        if any(tgt < 0)
            bad = enumIdx(find(tgt < 0, 1));
            error('ssit:pgf:negativeState', ...
                'Reaction %d (propensity ''%s'') can fire when species ''%s'' is %d, making it negative.', ...
                r, props{r}, species{bad}, e(enumIdx == bad));
        end
        if ~ismember(tgt', states', 'rows')
            states(:,end+1) = tgt; %#ok<AGROW>
            if size(states,2) > maxStates
                error('ssit:pgf:enumeratedUnbounded', ...
                    ['The enumerated species {%s} reach more than %d states. Choose ', ...
                     'pgfOptions.enumeratedSpecies explicitly or raise pgfOptions.maxEnumeratedStates.'], ...
                    strjoin(species(enumIdx), ', '), maxStates);
            end
        end
    end
    head = head + 1;
end
end

function B = randomStates(d, n, maxVal, enumIdx, enumStates, rs)
% n random nonnegative states for every enumerated state (d x n*nE).
nE = size(enumStates,2);
B = randi(rs, [0 maxVal], d, n*nE);
if ~isempty(enumIdx)
    B(enumIdx,:) = kron(enumStates, ones(1,n));
end
end

function opts = mergeOptions(src, opts)
defaults = ssit.pgf.defaultOptions();
if isempty(opts)
    if isPropOrField(src, 'pgfOptions') && ~isempty(src.pgfOptions)
        opts = src.pgfOptions;
    else
        opts = defaults;
    end
end
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(opts, names{i})
        opts.(names{i}) = defaults.(names{i});
    end
end
end

function tf = isPropOrField(obj, name)
tf = (isstruct(obj) && isfield(obj, name)) || (isobject(obj) && isprop(obj, name));
end

function p = probeParameters(parValues, rs)
p = parValues(:);
p(p == 0) = 1;
p = p .* (0.5 + rand(rs, size(p)));
end

function v = evalW(w, X, t, k, r, expr)
v = w(X, t, k);
v = v + zeros(1, size(X,2));    % constant propensities return a scalar
if ~isreal(v)
    % e.g. (negative)^3.1 * (t > t0) evaluates to a complex zero before t0
    if any(abs(imag(v)) > 1e-12 * max(1, abs(real(v))))
        error('ssit:pgf:nonPolynomial', ...
            'Reaction %d (propensity ''%s'') takes complex values at t = %g.', r, expr, t);
    end
    v = real(v);
end
if ~all(isfinite(v))
    error('ssit:pgf:nonPolynomial', ...
        'Reaction %d (propensity ''%s'') is not finite at some nonnegative integer states.', r, expr);
end
end

function X = latticePoints(K)
% All integer points 0 <= x <= K (numel(K) x prod(K+1)), first entry fastest.
d = numel(K);
n = prod(K + 1);
X = zeros(d, n);
idx = (0:n-1);
for i = 1:d
    X(i,:) = mod(idx, K(i) + 1);
    idx = floor(idx / (K(i) + 1));
end
end

function F = fallingFactorials(X, kOrders)
% F(q, n) = prod_i X(i,n)^(kOrders(q,i)) (falling factorial powers, orders <= 2).
F = ones(size(kOrders,1), size(X,2));
for q = 1:size(kOrders,1)
    for i = 1:size(X,1)
        switch kOrders(q,i)
            case 1
                F(q,:) = F(q,:) .* X(i,:);
            case 2
                F(q,:) = F(q,:) .* X(i,:) .* (X(i,:) - 1);
        end
    end
end
end

function w = makePropensityFunction(expr, species, parNames, inputs)
% Vectorised handle @(x,t,k): x is d x n (species), k is the parameter
% vector.  Inputs are substituted first (they may contain parameters).
e = char(expr);
for j = 1:size(inputs,1)
    e = regexprep(e, ['\<', inputs{j,1}, '\>'], ['(', inputs{j,2}, ')']);
end
% Whole-word substitution through placeholder tokens, so that a name can
% never match inside an already substituted one (e.g. species 'x' vs 'x(1,:)').
for i = 1:numel(species)
    e = regexprep(e, ['\<', species{i}, '\>'], sprintf('PGFSPECIESTOKEN%dEND', i));
end
for j = 1:numel(parNames)
    e = regexprep(e, ['\<', parNames{j}, '\>'], sprintf('PGFPARAMTOKEN%dEND', j));
end
e = regexprep(e, 'PGFSPECIESTOKEN(\d+)END', 'x($1,:)');
e = regexprep(e, 'PGFPARAMTOKEN(\d+)END', 'k($1)');
e = regexprep(e, '(?<!\.)([\*/\^])', '.$1');
w = str2func(['@(x,t,k)', e]);
end
