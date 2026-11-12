classdef pgfSolveTest < matlab.unittest.TestCase
    % Phase-2 tests for SSIT.solve(solver='pgf') (src/+ssit/+pgf).
    %
    % References: Julia fixtures (tests/test_data/pgf), analytic solutions
    % and SSIT's own FSP.  Expected accuracies come from a numpy/scipy mirror
    % of this code (operators vs Julia <= 2e-15; see the comments per test).
    %
    % Run from the tests folder:
    %   results = runtests('pgfSolveTest')

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

        function L = operatorAt(M, sg, t)
            % Generator matrix of the translated model on grid sg.
            model = ssit.pgf.translateModel(M, M.tSpan);
            ops = ssit.pgf.assembleOperators(model, sg);
            c = ssit.pgf.termCoefficients(ops, t, model.parameterValues);
            L = reshape(ops.Mflat * c, sg.N, sg.N);
        end

        function e = relErr(A, B)
            e = max(abs(A(:) - B(:))) / max(1, max(abs(B(:))));
        end
    end

    methods (Test)
        %% ---------------- Julia fixtures: operator and solution ----------------
        function birthDeathMatchesFixture(testCase)
            F = load(fullfile(testCase.FixtureDir, 'pgf_birthdeath.mat'));
            M = pgfSolveTest.makeModel({'x1'}, {'k';'g*x1'}, [1 -1], {'k',10;'g',1}, 0, F.times(:)');
            M.pgfOptions.mu = 4; M.pgfOptions.lowerBound = 0.8;
            M.pgfOptions.relTol = 1e-10; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            testCase.verifyLessThanOrEqual(max(abs(sol.grid.z - F.grid), [], 'all'), 1e-14);
            % Mirror: L vs Julia 5e-16, G vs fixture 6e-10, means vs analytic 2e-9.
            testCase.verifyLessThanOrEqual(pgfSolveTest.relErr(pgfSolveTest.operatorAt(M, sol.grid, 0), F.L), 1e-12);
            testCase.verifyLessThanOrEqual(max(abs(sol.G - F.G), [], 'all'), 1e-7);
            testCase.verifyLessThanOrEqual(max(abs(sol.means - F.means_exact(:)')), 1e-6);
        end

        function chainMatchesFixture(testCase)
            F = load(fullfile(testCase.FixtureDir, 'pgf_chain2.mat'));
            M = pgfSolveTest.makeModel({'x1','x2'}, {'b';'c*x1';'g*x1';'g*x2'}, ...
                [1 -1 -1 0; 0 1 0 -1], {'b',8;'c',2;'g',1}, [0;0], F.times(:)');
            M.pgfOptions.mu = [3 3]; M.pgfOptions.lowerBound = 0.85;
            M.pgfOptions.relTol = 1e-10; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            % Mirror: L vs Julia (Edward's rhs_chain_pgf!) 3e-16, G vs fixture 8e-10.
            testCase.verifyLessThanOrEqual(pgfSolveTest.relErr(pgfSolveTest.operatorAt(M, sol.grid, 0), F.L), 1e-12);
            testCase.verifyLessThanOrEqual(max(abs(sol.G - F.G), [], 'all'), 1e-7);
            % Discretisation error of mu = [3 3] vs the analytic means is 6e-4.
            testCase.verifyLessThanOrEqual(max(abs(sol.means - F.means_exact), [], 'all'), 1e-3);
        end

        function annihilationMatchesFixtureAndFsp(testCase)
            F = load(fullfile(testCase.FixtureDir, 'pgf_annihilation.mat'));
            M = pgfSolveTest.makeModel({'x1'}, {'k';'c*x1*(x1-1)'}, [1 -2], {'k',10;'c',0.1}, 0, F.times(:)');
            M.pgfOptions.mu = 5; M.pgfOptions.lowerBound = -1;
            M.pgfOptions.relTol = 1e-10; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            % Mirror: L vs Julia 2e-15, G vs fixture 9e-11, means vs FSP 2e-9.
            testCase.verifyLessThanOrEqual(pgfSolveTest.relErr(pgfSolveTest.operatorAt(M, sol.grid, 0), F.L), 1e-12);
            testCase.verifyLessThanOrEqual(max(abs(sol.G - F.G), [], 'all'), 1e-7);
            testCase.verifyLessThanOrEqual(max(abs(sol.means - F.means_fsp(:)')), 1e-6);
        end

        %% ---------------- analytic: nonzero initial state, auto domain ----------------
        function birthDeathFromNonzeroStateAuto(testCase)
            % x0 = 5: G = (1 - e + e z)^5 exp(m (z - 1)), e = exp(-g t), m = k/g (1 - e).
            k = 10; g = 1; T = linspace(0,5,11);
            M = pgfSolveTest.makeModel({'x1'}, {'k';'g*x1'}, [1 -1], {'k',k;'g',g}, 5, T);
            M.pgfOptions.mu = 4;     % lowerBound = 'auto' (default)
            M.pgfOptions.relTol = 1e-10; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            e = exp(-g*T); m = k/g*(1 - e);
            z = sol.grid.z;
            Gex = (1 - e + e.*z).^5 .* exp(m.*(z - 1));
            % Mirror: auto lowerBound 0.7993, G 3e-10, mean 3e-9, variance 9e-8.
            testCase.verifyEqual(sol.lowerBound, 0.7993, 'AbsTol', 1e-3);
            testCase.verifyLessThanOrEqual(max(abs(sol.G - Gex), [], 'all'), 1e-7);
            testCase.verifyLessThanOrEqual(max(abs(sol.means - (5*e + m))), 1e-6);
            testCase.verifyLessThanOrEqual(max(abs(squeeze(sol.covs)' - (5*e.*(1-e) + m))), 1e-5);
        end

        %% ---------------- FSP comparisons ----------------
        function timeVaryingInputMatchesFsp(testCase)
            T = linspace(0,6,13);
            M = pgfSolveTest.makeModel({'x1','x2'}, {'kr*I';'kt*x1';'g1*x1';'g2*x2'}, ...
                [1 -1 -1 0; 0 1 0 -1], {'kr',10;'w',1;'kt',2;'g1',1;'g2',0.5}, [0;0], T);
            M.inputExpressions = {'I','1+cos(w*t)'};
            M.fspOptions.fspTol = 1e-7;
            M.fspOptions.fspIntegratorRelTol = 1e-8;
            M.fspOptions.fspIntegratorAbsTol = 1e-12;
            M = M.solve(solver='fsp');
            M.pgfOptions.mu = 4; M.pgfOptions.lowerBound = 0.85;
            M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            Gfsp = ssit.pgf.fspToPgf(M.Solutions.fsp, sol.grid.z);
            meansFsp = zeros(2, numel(T));
            for it = 1:numel(T)
                X = M.Solutions.fsp{it}.p.data.subs - 1;
                meansFsp(:,it) = X' * M.Solutions.fsp{it}.p.data.vals;
            end
            % Mirror (scipy FSP, 3e-15 lost mass): G 1.7e-6, means 9e-8.
            testCase.verifyLessThanOrEqual(max(abs(sol.G - Gfsp), [], 'all'), 1e-4);
            testCase.verifyLessThanOrEqual(max(abs(sol.means - meansFsp), [], 'all'), 1e-3);
        end

        function bimolecularMatchesFsp(testCase)
            % Mixed bimolecular A + B -> 0 needs a fine grid: the mirror's
            % mean error vs FSP is 0.6 / 0.1 / 6e-3 / 1e-5 / 2e-8 for
            % mu = 3..7 on [0.5, 1]^2.
            T = linspace(0,4,5);
            M = pgfSolveTest.makeModel({'A','B'}, {'kA';'kB';'c*A*B';'g*A';'g*B'}, ...
                [1 0 -1 -1 0; 0 1 -1 0 -1], {'kA',10;'kB',8;'c',0.1;'g',0.5}, [0;0], T);
            M.fspOptions.fspTol = 1e-7;
            M.fspOptions.fspIntegratorRelTol = 1e-8;
            M.fspOptions.fspIntegratorAbsTol = 1e-12;
            M = M.solve(solver='fsp');
            M.pgfOptions.mu = 6; M.pgfOptions.lowerBound = 0.5;
            M.pgfOptions.relTol = 1e-9; M.pgfOptions.absTol = 1e-12;
            M = M.solve(solver='pgf');
            sol = M.Solutions.pgf;
            Gfsp = ssit.pgf.fspToPgf(M.Solutions.fsp, sol.grid.z);
            testCase.verifyLessThanOrEqual(max(abs(sol.G - Gfsp), [], 'all'), 1e-4);
        end

        %% ---------------- error paths ----------------
        function nonPolynomialPropensityErrors(testCase)
            M = pgfSolveTest.makeModel({'x1'}, {'k';'g*x1/(1+x1)'}, [1 -1], {'k',10;'g',1}, 0, linspace(0,5,6));
            testCase.verifyError(@() M.solve(solver='pgf'), 'ssit:pgf:nonPolynomial');
            M.propensityFunctions = {'k';'g*x1^3'};
            testCase.verifyError(@() M.solve(solver='pgf'), 'ssit:pgf:nonPolynomial');
        end

        function negativeStateErrors(testCase)
            M = pgfSolveTest.makeModel({'x1'}, {'k';'g'}, [1 -1], {'k',10;'g',1}, 0, linspace(0,5,6));
            testCase.verifyError(@() M.solve(solver='pgf'), 'ssit:pgf:negativeState');
        end

        function secondOrderNeedsFullDomain(testCase)
            M = pgfSolveTest.makeModel({'x1'}, {'k';'c*x1*(x1-1)'}, [1 -2], {'k',10;'c',0.1}, 0, linspace(0,5,6));
            M.pgfOptions.lowerBound = 0.8;
            testCase.verifyError(@() M.solve(solver='pgf'), 'ssit:pgf:secondOrderDomain');
            M.pgfOptions.lowerBound = 'auto';
            M = M.solve(solver='pgf');
            testCase.verifyEqual(M.Solutions.pgf.lowerBound, -1);
        end

        function unstableOperatorErrors(testCase)
            % 0 -> A, 0 -> B, A + B -> 0 on [0.9,1]^2 with mu = 5: max Re(eig) = +2.7.
            M = pgfSolveTest.makeModel({'A','B'}, {'kA';'kB';'c*A*B'}, [1 0 -1; 0 1 -1], ...
                {'kA',10;'kB',8;'c',0.1}, [0;0], linspace(0,2,5));
            M.pgfOptions.mu = 5; M.pgfOptions.lowerBound = 0.9;
            testCase.verifyError(@() M.solve(solver='pgf'), 'ssit:pgf:unstable');
        end

        function zeroParameterKeepsStructure(testCase)
            % The operator structure must not depend on current values
            % (a fit may start from g = 0).
            M = pgfSolveTest.makeModel({'x1'}, {'k';'g*x1'}, [1 -1], {'k',10;'g',0}, 0, linspace(0,5,6));
            model = ssit.pgf.translateModel(M, M.tSpan);
            testCase.verifyEqual(model.rxn(2).k, 1);
        end
    end
end
