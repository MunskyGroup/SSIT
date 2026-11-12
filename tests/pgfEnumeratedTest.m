classdef pgfEnumeratedTest < matlab.unittest.TestCase
    % Phase-3 tests for the PGF solver with enumerated (gene-type) species
    % and steady-state initial conditions (src/+ssit/+pgf).
    %
    % References: Edward's GR-DUSP1 fixture (tests/test_data/pgf), analytic
    % telegraph-model steady state, and SSIT's own FSP.  Expected accuracies
    % come from a numpy/scipy mirror of this code (see comments per test).
    %
    % Run from the tests folder:
    %   results = runtests('pgfEnumeratedTest')

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
        function M = makeModel(species, props, S, pars, x0, tSpan)
            M = SSIT;
            M.species = species;
            M.initialCondition = x0;
            M.propensityFunctions = props;
            M.stoichiometry = S;
            M.parameters = pars;
            M.tSpan = tSpan;
        end

        function M = burstingGene(tSpan)
            M = pgfEnumeratedTest.makeModel({'geneState';'rna'}, ...
                {'kon*(1-geneState)';'koff*geneState';'kr*geneState';'gr*rna'}, ...
                [1 -1 0 0; 0 0 1 -1], {'kon',1;'koff',1;'kr',10;'gr',1}, [0;0], tSpan);
        end

        function M = grdusp1(pars, tSpan)
            M = pgfEnumeratedTest.makeModel({'A';'B';'C';'D';'E';'F'}, ...
                {'kg';'gc*A';'kcn*A';'knc*B';'gn*B';'kon*B*C';'koff*D';'kr*D';'kex*E';'kdeg*F'}, ...
                [1 -1 -1  1  0  0  0 0  0  0;
                 0  0  1 -1 -1 -1  1 0  0  0;
                 0  0  0  0  0 -1  1 0  0  0;
                 0  0  0  0  0  1 -1 0  0  0;
                 0  0  0  0  0  0  0 1 -1  0;
                 0  0  0  0  0  0  0 0  1 -1], ...
                pars, [0;0;1;0;0;0], tSpan);
            M.inputExpressions = {'kcn','kcn_pre+(kcn_post-kcn_pre)*(t>0)'};
        end

        function [means, covs] = fspMoments(fsp)
            Nt = numel(fsp);
            d = size(fsp{1}.p.data.subs, 2);
            means = zeros(d, Nt); covs = zeros(d, d, Nt);
            for it = 1:Nt
                X = fsp{it}.p.data.subs - 1;
                p = fsp{it}.p.data.vals;
                means(:,it) = X' * p;
                covs(:,:,it) = X' * (p .* X) - means(:,it) * means(:,it)';
            end
        end

        function M = solveFsp(M)
            M.fspOptions.fspTol = 1e-8;
            M.fspOptions.fspIntegratorRelTol = 1e-9;
            M.fspOptions.fspIntegratorAbsTol = 1e-13;
            M = M.solve(solver='fsp');
        end

        function verifyAgainstFsp(testCase, M, tolG, tolMean, tolCovRel)
            sol = M.Solutions.pgf;
            Y = find(ismember(M.species, sol.pgfSpecies));
            Gfsp = ssit.pgf.fspToPgf(M.Solutions.fsp, sol.grid.z, Y);
            [mF, cF] = pgfEnumeratedTest.fspMoments(M.Solutions.fsp);
            testCase.verifyLessThanOrEqual(max(abs(sol.Gtot - Gfsp), [], 'all'), tolG);
            testCase.verifyLessThanOrEqual(max(abs(sol.means - mF), [], 'all'), tolMean);
            testCase.verifyLessThanOrEqual(max(abs(sol.covs - cF), [], 'all') / max(abs(cF), [], 'all'), tolCovRel);
        end
    end

    methods (Test)
        %% ---------------- automatic detection ----------------
        function autoDetectionOnExampleModels(testCase)
            T = linspace(0,10,5);
            centralDogma = pgfEnumeratedTest.makeModel({'rna';'protein'}, {'kr';'gr*rna';'kp*rna';'gp*protein'}, ...
                [1 -1 0 0; 0 0 1 -1], {'kr',10;'gr',1;'kp',1;'gp',0.1}, [0;0], T);
            m = ssit.pgf.translateModel(centralDogma, T);
            testCase.verifyEmpty(m.enumIdx);

            m = ssit.pgf.translateModel(pgfEnumeratedTest.burstingGene(T), T);   % bounded by its propensity
            testCase.verifyEqual(m.species(m.enumIdx), {'geneState'});
            testCase.verifyEqual(m.enumStates, [0 1]);

            offOn = pgfEnumeratedTest.makeModel({'offGene';'onGene';'mRNA'}, ...
                {'kon*offGene';'koff*onGene';'kr*onGene';'gr*mRNA'}, [-1 1 0 0; 1 -1 0 0; 0 0 1 -1], ...
                {'kon',0.5;'koff',1;'kr',20;'gr',1}, [1;0;0], T);
            m = ssit.pgf.translateModel(offOn, T);                               % conservation law
            testCase.verifyEqual(m.species(m.enumIdx), {'offGene','onGene'});
            testCase.verifyEqual(m.enumStates, [1 0; 0 1]);

            pars = [{'kg';'gc';'kcn_pre';'kcn_post';'knc';'gn';'kon';'koff';'kr';'kex';'kdeg'}, num2cell(ones(11,1))];
            m = ssit.pgf.translateModel(pgfEnumeratedTest.grdusp1(pars, T), T);
            testCase.verifyEqual(m.species(m.enumIdx), {'C','D'});
            testCase.verifyTrue(m.timeDependent);
            testCase.verifyFalse(any(m.secondOrder | m.mixedOrder));   % kon*B*C is first order in B

            stl1 = SSIT('STL1_4state');
            m = ssit.pgf.translateModel(stl1, linspace(0,30,5));
            testCase.verifyEqual(m.species(m.enumIdx), {'g1','g2','g3','g4'});
            testCase.verifyEqual(size(m.enumStates,2), 4);
        end

        %% ---------------- Edward's GR-DUSP1 model from SSIT strings ----------------
        function grdusp1MatchesFixture(testCase)
            F = load(fullfile(testCase.FixtureDir, 'pgf_grdusp1.mat'));
            names = {'kg';'gc';'kcn_pre';'kcn_post';'knc';'gn';'kon';'koff';'kr';'kex';'kdeg'};
            pars = [names, num2cell(F.truth.parameters(:))];
            tSpan = [0, F.post_times(:)' - F.t_switch];
            M = pgfEnumeratedTest.grdusp1(pars, tSpan);
            M.fspOptions.initApproxSS = true;           % steady state at t = 0 (pre-Dex)
            M.pgfOptions.mu = [4 2 2 2];                 % per PGF species A, B, E, F
            M.pgfOptions.lowerBound = 0.95;
            M.pgfOptions.relTol = 1e-10; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            N = sol.grid.N;
            testCase.verifyEqual(sol.enumeratedSpecies, {'C','D'});
            testCase.verifyLessThanOrEqual(max(abs(sol.grid.z - F.grid), [], 'all'), 1e-14);
            % Operators (Edward's hand-built blocks) -- mirror: 2e-15.
            model = ssit.pgf.translateModel(M, tSpan);
            ops = ssit.pgf.assembleOperators(model, sol.grid);
            Lpre = ssit.pgf.assembleGenerator(ops, 0, model.parameterValues);
            Lpost = ssit.pgf.assembleGenerator(ops, 1, model.parameterValues);
            rel = @(a,b) max(abs(a(:)-b(:))) / max(1, max(abs(b(:))));
            testCase.verifyLessThanOrEqual(rel(Lpre * F.v_test(:), F.truth.L_pre_times_v(:)), 1e-12);
            testCase.verifyLessThanOrEqual(rel(Lpost * F.v_test(:), F.truth.L_post_times_v(:)), 1e-12);
            % Steady state (mirror 7e-12) and the 9 snapshots (mirror 3.5e-10).
            testCase.verifyLessThanOrEqual(max(abs(sol.G(:,1) - F.truth.G_steady(:))), 1e-9);
            testCase.verifyLessThanOrEqual(max(abs(sol.Gtot - F.truth.gtot), [], 'all'), 1e-7);
            testCase.verifyEqual(size(sol.G, 1), 2*N);
        end

        %% ---------------- general models vs FSP ----------------
        function burstingGeneMatchesFsp(testCase)
            % Mirror: Gtot 1.3e-9, means 1.1e-8, covs 6.7e-8 (max cov 13).
            M = pgfEnumeratedTest.burstingGene(linspace(0,10,11));
            M = pgfEnumeratedTest.solveFsp(M);
            M.pgfOptions.mu = 4; M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            testCase.verifyEqual(M.Solutions.pgf.enumeratedSpecies, {'geneState'});
            pgfEnumeratedTest.verifyAgainstFsp(testCase, M, 1e-6, 1e-5, 1e-5);
        end

        function burstingGeneWithoutEnumeration(testCase)
            % enumeratedSpecies = {}: geneState becomes a PGF variable on
            % [-1,1]; still correct but much less accurate (mirror: means
            % 1.2e-4, covs 1.3e-2 absolute) -- the reason to enumerate.
            M = pgfEnumeratedTest.burstingGene(linspace(0,10,11));
            M = pgfEnumeratedTest.solveFsp(M);
            M.pgfOptions.enumeratedSpecies = {};
            M.pgfOptions.mu = [2 4]; M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            testCase.verifyEmpty(M.Solutions.pgf.enumeratedSpecies);
            pgfEnumeratedTest.verifyAgainstFsp(testCase, M, 1e-5, 1e-3, 1e-2);
        end

        function threeStateGeneMatchesFsp(testCase)
            % example_01 model: offGene + onGene = 1.  Mirror: Gtot 2.2e-9, means 3.7e-8.
            M = pgfEnumeratedTest.makeModel({'offGene';'onGene';'mRNA'}, ...
                {'kon*offGene';'koff*onGene';'kr*onGene';'gr*mRNA'}, [-1 1 0 0; 1 -1 0 0; 0 0 1 -1], ...
                {'kon',0.5;'koff',1;'kr',20;'gr',1}, [1;0;0], linspace(0,8,9));
            M = pgfEnumeratedTest.solveFsp(M);
            M.pgfOptions.mu = 4; M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            pgfEnumeratedTest.verifyAgainstFsp(testCase, M, 1e-6, 1e-5, 1e-5);
        end

        function stl1PresetMatchesFsp(testCase)
            % SSIT preset STL1_4state: 4 enumerated gene states, Hog1 input
            % with max() and a complex-zero branch before t0.  Mirror (mu = 5):
            % Gtot 1.3e-9, means 5.9e-8, covs 2.8e-5 (max cov 262).
            M = SSIT('STL1_4state');
            M.tSpan = linspace(0,30,16);
            M = pgfEnumeratedTest.solveFsp(M);
            M.pgfOptions.mu = 5; M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            testCase.verifyEqual(M.Solutions.pgf.enumeratedSpecies, {'g1','g2','g3','g4'});
            pgfEnumeratedTest.verifyAgainstFsp(testCase, M, 1e-5, 1e-3, 1e-4);
        end

        %% ---------------- steady state ----------------
        function telegraphSteadyStateAnalytic(testCase)
            % Telegraph model: mean = kr kon/((kon+koff) gr),
            % Fano = 1 + kr koff/((kon+koff)(kon+koff+gr)), P(on) = kon/(kon+koff).
            % Mirror: exact to 10 digits; stationary to 2e-15.
            M = pgfEnumeratedTest.burstingGene([0 1]);
            M.fspOptions.initApproxSS = true;
            M.pgfOptions.mu = 4; M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            kon = 1; koff = 1; kr = 10; gr = 1;
            testCase.verifyEqual(sol.means(2,1), kr*kon/((kon+koff)*gr), 'AbsTol', 1e-8);
            testCase.verifyEqual(sol.covs(2,2,1)/sol.means(2,1), 1 + kr*koff/((kon+koff)*(kon+koff+gr)), 'AbsTol', 1e-8);
            testCase.verifyEqual(sol.enumeratedProbabilities(2,1), kon/(kon+koff), 'AbsTol', 1e-10);
            testCase.verifyEqual(sol.means(1,1), kon/(kon+koff), 'AbsTol', 1e-10);   % gene mean = P(on)
            testCase.verifyLessThanOrEqual(max(abs(sol.G(:,2) - sol.G(:,1))), 1e-8);
        end

        %% ---------------- reuse of the setup between solves ----------------
        function setupReusedBetweenSolves(testCase)
            M = pgfEnumeratedTest.burstingGene(linspace(0,5,6));
            M.pgfOptions.mu = 4; M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            c1 = M.pgfCache;
            testCase.verifyNotEmpty(c1);
            G1 = M.Solutions.pgf.G;

            % Same model: same grid and domain, same solution.
            M = M.solve(solver='pgf');
            testCase.verifyEqual(M.pgfCache.lower, c1.lower);
            testCase.verifyEqual(M.Solutions.pgf.G, G1, 'AbsTol', 1e-12);

            % New parameter value (mean 5 -> 6, domain still fits): cached
            % setup, and the same moments as a solve from scratch.
            M.parameters{3,2} = 12;
            M = M.solve(solver='pgf');
            testCase.verifyEqual(M.pgfCache.lower, c1.lower);
            Mfresh = M;
            Mfresh.pgfCache = [];
            Mfresh = Mfresh.solve(solver='pgf');
            testCase.verifyEqual(M.Solutions.pgf.means, Mfresh.Solutions.pgf.means, 'AbsTol', 1e-5);

            % Changed structure: the setup is rebuilt.
            M.propensityFunctions{4} = 'gr*rna*(1+0*geneState)';
            M = M.solve(solver='pgf');
            testCase.verifyNotEqual(M.pgfCache.key, c1.key);
        end

        %% ---------------- errors ----------------
        function unboundedEnumerationErrors(testCase)
            M = pgfEnumeratedTest.burstingGene(linspace(0,5,6));
            M.pgfOptions.enumeratedSpecies = {'rna'};
            M.pgfOptions.maxEnumeratedStates = 20;
            testCase.verifyError(@() M.solve(solver='pgf'), 'ssit:pgf:enumeratedUnbounded');
        end

        function allEnumeratedErrors(testCase)
            M = pgfEnumeratedTest.makeModel({'offGene';'onGene'}, {'kon*offGene';'koff*onGene'}, ...
                [-1 1; 1 -1], {'kon',1;'koff',2}, [1;0], linspace(0,5,6));
            testCase.verifyError(@() M.solve(solver='pgf'), 'ssit:pgf:noPgfSpecies');
        end

        function hillFunctionOfEnumeratedSpeciesIsAllowed(testCase)
            % Propensities only need to be polynomial in the PGF species: a
            % Hill function of an enumerated (bounded) species is fine.
            M = pgfEnumeratedTest.makeModel({'g';'rna'}, ...
                {'kon*(1-g)';'koff*g';'kr*(0.1+g^2/(0.5+g^2))';'gr*rna'}, ...
                [1 -1 0 0; 0 0 1 -1], {'kon',1;'koff',2;'kr',15;'gr',1}, [0;0], linspace(0,6,7));
            M = pgfEnumeratedTest.solveFsp(M);
            M.pgfOptions.mu = 4; M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            pgfEnumeratedTest.verifyAgainstFsp(testCase, M, 1e-6, 1e-5, 1e-5);
        end
    end
end
