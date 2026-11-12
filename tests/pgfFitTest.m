classdef pgfFitTest < matlab.unittest.TestCase
    % Phase-4 tests: fitting CME parameters to data with the PGF solver
    % (ssit.pgf.fitSetup / residuals / logLikelihood / levenbergMarquardt
    % and the solutionScheme = 'pgf' branches of SSIT.computeLikelihood and
    % SSIT.maximizeLikelihood).
    %
    % References: Edward's GR-DUSP1 objective (tests/test_data/pgf), exact
    % Poisson maximum-likelihood estimates for birth-death data
    % (bd_*_counts.csv, written by the numpy mirror with seed 20261109), and
    % SSIT's FSP for the grouping of enumerated / unobserved species.
    % Expected accuracies come from a numpy/scipy mirror (comments per test).
    %
    % Run from the tests folder:
    %   results = runtests('pgfFitTest')

    properties
        FixtureDir
    end

    methods (TestClassSetup)
        function setPaths(testCase)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fullfile(here, '..', 'src')));
            testCase.FixtureDir = fullfile(here, 'test_data', 'pgf');
        end
    end

    methods (Static)
        function M = makeModel(species, props, S, pars, x0)
            M = SSIT;
            M.species = species;
            M.initialCondition = x0;
            M.propensityFunctions = props;
            M.stoichiometry = S;
            M.parameters = pars;
        end

        function M = grdusp1(dataFile, pars)
            % Edward's GR-DUSP1 model in hours, Dex added at t = 150 h; the
            % pre-Dex steady state is the initial condition at t = 150.
            M = pgfFitTest.makeModel({'A';'B';'C';'D';'E';'F'}, ...
                {'kg';'gc*A';'kcn*A';'knc*B';'gn*B';'kon*B*C';'koff*D';'kr*D';'kex*E';'kdeg*F'}, ...
                [1 -1 -1  1  0  0  0 0  0  0;
                 0  0  1 -1 -1 -1  1 0  0  0;
                 0  0  0  0  0 -1  1 0  0  0;
                 0  0  0  0  0  1 -1 0  0  0;
                 0  0  0  0  0  0  0 1 -1  0;
                 0  0  0  0  0  0  0 0  1 -1], pars, [0;0;1;0;0;0]);
            M.inputExpressions = {'kcn','kcn_pre+(kcn_post-kcn_pre)*(t>150)'};
            M.initialTime = 150;
            M.fspOptions.initApproxSS = true;
            % Only A, B, E, F are linked (as in Edward's fit); C and D are
            % in the file but unobserved.  The 100 h snapshot precedes the
            % initial time and is dropped.
            M = M.loadData(dataFile, {'A','A';'B','B';'E','E';'F','F'}, {'time',100,'>'});
            M.solutionScheme = 'pgf';
            M.pgfOptions.mu = [4 2 2 2];
            M.pgfOptions.lowerBound = 0.95;
        end

        function M = birthDeath(dataFile, k0)
            M = pgfFitTest.makeModel({'x1'}, {'k';'g*x1'}, [1 -1], {'k',k0;'g',1}, 0);
            M = M.loadData(dataFile, {'x1','x1'});
            M.solutionScheme = 'pgf';
        end

        function [pHat, sd] = poissonMle(times, counts)
            % Exact MLE of (k, g) for x(t) ~ Poisson(k/g (1 - exp(-g t)))
            % from x(0) = 0, and the SD of log(k), log(g) (inverse Hessian).
            lam = @(p, t) exp(p(1)) / exp(p(2)) * (1 - exp(-exp(p(2)) * t));
            tt = unique(times);
            S = arrayfun(@(t) sum(counts(times == t)), tt);
            n = arrayfun(@(t) sum(times == t), tt);
            nll = @(p) -sum(S .* log(lam(p, tt)) - n .* lam(p, tt));
            p = fminsearch(nll, log([5; 1]), optimset('TolX', 1e-12, 'TolFun', 1e-12, 'MaxFunEvals', 1e4, 'MaxIter', 1e4));
            h = 1e-4; H = zeros(2);
            for i = 1:2
                for j = 1:2
                    ei = zeros(2,1); ei(i) = h; ej = zeros(2,1); ej(j) = h;
                    H(i,j) = (nll(p+ei+ej) - nll(p+ei-ej) - nll(p-ei+ej) + nll(p-ei-ej)) / (4*h*h);
                end
            end
            pHat = exp(p);
            sd = sqrt(diag(inv(H)));
        end
    end

    methods (Test)
        %% ---------------- Edward's objective ----------------
        function grdusp1MseMatchesEdward(testCase)
            % Mirror: data PGF 6e-14, -logL(truth) = 1.03974e-5,
            % -logL(ones) = 0.224404 (fixture 1.0398e-5 and 0.22440).
            F = load(fullfile(testCase.FixtureDir, 'pgf_grdusp1.mat'));
            names = {'kg';'gc';'kcn_pre';'kcn_post';'knc';'gn';'kon';'koff';'kr';'kex';'kdeg'};
            M = pgfFitTest.grdusp1(fullfile(testCase.FixtureDir, 'grdusp1_steady_counts.csv'), ...
                [names, num2cell(F.truth.parameters(:))]);
            M.pgfOptions.loss = 'mse';
            fp = ssit.pgf.fitSetup(M);
            testCase.verifyEqual(fp.observedSpecies, {'A','B','E','F'});
            testCase.verifyEqual(fp.nG, 1);
            testCase.verifyEqual(fp.fitTimes, F.snapshot_times(:)', 'AbsTol', 1e-9);
            testCase.verifyLessThanOrEqual(max(abs(fp.prob.sg.z - F.grid), [], 'all'), 1e-14);
            for k = 1:9
                testCase.verifyLessThanOrEqual(max(abs(fp.data(k).m - F.ref_pgf(:, k+1))), 1e-12);
            end
            testCase.verifyEqual(-ssit.pgf.logLikelihood(fp), F.truth.mse_objective, 'RelTol', 1e-3);
            testCase.verifyEqual(-ssit.pgf.logLikelihood(fp, ones(11,1)), F.ones.mse_objective, 'RelTol', 1e-3);
            % SSIT wrapper gives the same value.
            testCase.verifyEqual(M.computeLikelihood(), ssit.pgf.logLikelihood(fp), 'RelTol', 1e-10);
        end

        function grdusp1GlsErrorFloor(testCase)
            % 'gls' with errorFloor 'auto': the PGF on [4 2 2 2] differs from
            % [5 3 3 3] by 3.75e-4 RMS at truth (mirror), the error floor.
            % -logL(truth) = 136.6, -logL(ones) = 8.05e6 (mirror).
            F = load(fullfile(testCase.FixtureDir, 'pgf_grdusp1.mat'));
            names = {'kg';'gc';'kcn_pre';'kcn_post';'knc';'gn';'kon';'koff';'kr';'kex';'kdeg'};
            M = pgfFitTest.grdusp1(fullfile(testCase.FixtureDir, 'grdusp1_steady_counts.csv'), ...
                [names, num2cell(F.truth.parameters(:))]);
            fp = ssit.pgf.fitSetup(M);
            testCase.verifyEqual(fp.loss, 'gls');
            testCase.verifyEqual(fp.errorFloor, 3.75e-4, 'RelTol', 0.02);
            testCase.verifyEqual(-ssit.pgf.logLikelihood(fp), 136.6, 'RelTol', 0.02);
            testCase.verifyEqual(-ssit.pgf.logLikelihood(fp, ones(11,1)), 8.05e6, 'RelTol', 0.02);
        end

        %% ---------------- exact maximum likelihood ----------------
        function birthDeathSteadyStateMatchesMle(testCase)
            % 2000 cells from Poisson(10) at steady state; fit k (g = 1).
            % The MLE is the sample mean (9.9825, SE 0.0706).  Mirror: gls
            % 9.97744 (0.07 SE away), SD from J'J 0.0693.
            file = fullfile(testCase.FixtureDir, 'bd_steady_counts.csv');
            D = readmatrix(file);
            mle = mean(D(:,2));
            se = sqrt(mle / size(D,1));
            M = pgfFitTest.birthDeath(file, 5);
            M.fspOptions.initApproxSS = true;
            M.fittingOptions.modelVarsToFit = 1;
            opts = struct('Display', 'off');
            [k, ~, info] = M.maximizeLikelihood([], opts, 'levenbergMarquardt', returnType = 'parameters');
            testCase.verifyLessThanOrEqual(abs(k - mle), 0.25 * se);
            testCase.verifyEqual(k * sqrt(1 / info.fimLog), se, 'RelTol', 0.1);
            % fminsearch through the same objective lands on the same point.
            k2 = M.maximizeLikelihood([], optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-12), ...
                'fminsearch', returnType = 'parameters');
            testCase.verifyEqual(k2, k, 'RelTol', 1e-4);
            % returnType 'ssit' stores the estimate in the model.
            M2 = M.maximizeLikelihood([], opts, 'levenbergMarquardt');
            testCase.verifyEqual(M2.parameters{1,2}, k, 'RelTol', 1e-12);
        end

        function birthDeathTimeCourseMatchesMle(testCase)
            % 4 x 500 cells at t = 1, 2, 4, 8 from x(0) = 0 (k = 10, g = 0.5).
            % Exact MLE k, g = 9.6398, 0.47692 (SD of logs 0.0162, 0.0238).
            % Mirror |log(est/MLE)|/SD: gls 0.04 / 0.09, weighted 0.22 /
            % 0.12, mse 0.48 / 0.37; SD from J'J 0.0155 / 0.0230.
            file = fullfile(testCase.FixtureDir, 'bd_timecourse_counts.csv');
            D = readmatrix(file);
            [mle, sd] = pgfFitTest.poissonMle(D(:,1), D(:,2));
            opts = struct('Display', 'off');
            tol = struct('gls', 0.25, 'weighted', 0.5, 'mse', 1);
            for loss = {'gls', 'weighted', 'mse'}
                M = pgfFitTest.birthDeath(file, 5);
                M.pgfOptions.loss = loss{1};
                [p, ~, info] = M.maximizeLikelihood([], opts, 'levenbergMarquardt', returnType = 'parameters');
                testCase.verifyLessThanOrEqual(abs(log(p(:) ./ mle)) ./ sd, tol.(loss{1}), loss{1});
                if strcmp(loss{1}, 'gls')
                    testCase.verifyEqual(sqrt(diag(inv(info.fimLog))), sd, 'RelTol', 0.15);
                end
            end
        end

        function metropolisHastingsRuns(testCase)
            % Smoke test: MH samples the 'gls' approximate posterior.
            file = fullfile(testCase.FixtureDir, 'bd_steady_counts.csv');
            M = pgfFitTest.birthDeath(file, 10);
            M.fspOptions.initApproxSS = true;
            M.fittingOptions.modelVarsToFit = 1;
            opts = struct('numberOfSamples', 40, 'burnIn', 0, 'progress', false, ...
                'proposalDistribution', @(x) x + 0.007*randn(size(x)), ...
                'saveFile', [tempname, '.mat'], 'computeESS', false);
            [~, ~, res] = M.maximizeLikelihood([], opts, 'MetropolisHastings', returnType = 'parameters');
            testCase.verifyNumElements(res.mhSamples, 40);
            testCase.verifyTrue(all(isfinite(res.mhValue)));
            testCase.verifyLessThan(max(abs(exp(res.mhSamples) - 9.98)), 0.5);
        end

        %% ---------------- groups of enumerated states, unobserved species ----------------
        function groupsAndUnobservedSpeciesMatchFsp(testCase)
            % "Data" = the FSP distribution itself (tensor values = probabilities),
            % so model - data is the PGF discretisation error only (it falls
            % with mu: protein-only 9e-4 at mu = 4, 2.5e-5 at mu = 5).
            M = pgfFitTest.makeModel({'geneState';'rna';'protein'}, ...
                {'kon*(1-geneState)';'koff*geneState';'kr*geneState';'gr*rna';'kp*rna';'gp*protein'}, ...
                [1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1], ...
                {'kon',1;'koff',2;'kr',10;'gr',1;'kp',2;'gp',0.5}, [0;0;0]);
            T = [0 1 3 6];
            M.tSpan = T;
            M.fspOptions.fspTol = 1e-9;
            M.fspOptions.fspIntegratorRelTol = 1e-10;
            M.fspOptions.fspIntegratorAbsTol = 1e-14;
            M = M.solve(solver='fsp');
            subs = []; vals = [];
            for it = 1:numel(T)
                d = M.Solutions.fsp{it}.p.data;
                subs = [subs; it*ones(size(d.subs,1),1), d.subs]; %#ok<AGROW>
                vals = [vals; d.vals]; %#ok<AGROW>
            end
            tensor = sptensor(subs, vals, [numel(T), max(subs(:,2:end), [], 1)]);
            M.solutionScheme = 'pgf';
            M.pgfOptions.mu = 5;
            M.pgfOptions.lowerBound = [0.4 0.83];
            M.pgfOptions.loss = 'mse';
            M.fittingOptions.timesToFit = [2 3 4];
            % Mirror max |model - data|: 2.5e-5, 2.7e-7, 1.7e-5, 2.5e-5, 1.9e-7.
            cases = {{'geneState','rna','protein'}, 2, 1e-4;
                     {'rna'},                       1, 2e-6;
                     {'geneState','protein'},       2, 1e-4;
                     {'protein'},                   1, 1e-4;
                     {'geneState'},                 2, 1e-6};
            for c = 1:size(cases, 1)
                obs = cases{c,1};
                M.dataSet = struct('times', T, 'tensorSpecies', {M.species'}, ...
                    'linkedSpecies', dictionary(string(obs), string(obs)), ...
                    'app', struct('DataLoadingAndFittingTabOutputs', struct('dataTensor', tensor)));
                fp = ssit.pgf.fitSetup(M);
                testCase.verifyEqual(fp.nG, cases{c,2}, strjoin(obs, ','));
                r = ssit.pgf.residuals(fp) * sqrt(3 * fp.M * fp.nG);     % undo the MSE scaling
                testCase.verifyLessThanOrEqual(max(abs(r)), cases{c,3}, strjoin(obs, ','));
            end
        end

        %% ---------------- errors ----------------
        function errors(testCase)
            M = pgfFitTest.makeModel({'x1'}, {'k';'g*x1'}, [1 -1], {'k',10;'g',1}, 0);
            M.solutionScheme = 'pgf';
            testCase.verifyError(@() ssit.pgf.fitSetup(M), 'ssit:pgf:noData');
            M = M.loadData(fullfile(testCase.FixtureDir, 'bd_timecourse_counts.csv'), {'x1','x1'});
            M.initialTime = 3;                       % data at t = 1, 2 precede it
            testCase.verifyError(@() ssit.pgf.fitSetup(M), 'ssit:pgf:data');
            M.initialTime = 0;
            M.pgfOptions.loss = 'likelihood';
            testCase.verifyError(@() ssit.pgf.fitSetup(M), 'ssit:pgf:options');
        end
    end
end
