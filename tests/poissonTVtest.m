classdef poissonTVtest < matlab.unittest.TestCase
    properties
        TvPoiss
        TvPoissSolution
        PoissTVODE
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(testCase1)
            addpath('../CommandLine')
            %% Test Case - a simple Poisson model with time varying production
            testCase1.TvPoiss = SSIT;
            testCase1.TvPoiss.species = {'rna'};
            testCase1.TvPoiss.initialCondition = 0;
            testCase1.TvPoiss.stoichiometry = [1,-1];
            testCase1.TvPoiss.parameters = ({'kr',10;'gr',1});
            testCase1.TvPoiss.tSpan = linspace(0,2,21);
            testCase1.TvPoiss.fspOptions.fspTol = 1e-5;

            testCase1.TvPoiss.propensityFunctions = {'kr*Ig';'gr*rna'};
            testCase1.TvPoiss.inputExpressions = {'Ig','t>1'};
            testCase1.TvPoiss = testCase1.TvPoiss.formPropensitiesGeneral('PoissTV');

            [testCase1.TvPoissSolution,testCase1.TvPoiss.fspOptions.bounds] = testCase1.TvPoiss.solve;
            tic
            [testCase1.TvPoissSolution,testCase1.TvPoiss.fspOptions.bounds] = testCase1.TvPoiss.solve(testCase1.TvPoissSolution.stateSpace);
            testCase1.TvPoissSolution.time = toc;

            %% ODE model for TV Poisson process
            testCase1.PoissTVODE = testCase1.TvPoiss;
            testCase1.PoissTVODE.solutionScheme = 'ode';
            testCase1.PoissTVODE = testCase1.PoissTVODE.formPropensitiesGeneral('PoissTVODE');  

         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)

        function TimeVaryingPoissonMeansSuccess(testCase)
            % In this test, we check that the 1D Time Varying Poisson model
            % generates a solution with the correct mean versus time.
            t = testCase.TvPoiss.tSpan;
            tp = max(0,t-1);
            
            mn = testCase.TvPoiss.parameters{1,2}/testCase.TvPoiss.parameters{2,2}*...
                (1-exp(-testCase.TvPoiss.parameters{2,2}*tp));
            
            fspSoln = testCase.TvPoissSolution.fsp;
            fspMn = NaN*mn;
            for i = 1:length(fspSoln)
                n = size(testCase.TvPoissSolution.fsp{i}.p.data,1);
                fspMn(i) = [0:n-1]*double(testCase.TvPoissSolution.fsp{i}.p.data);
            end

            diff = abs(fspMn-mn);
            relDiff = sum(diff(mn>0)./mn(mn>0))/length(fspSoln);

            testCase.verifyEqual(relDiff<0.01, true, ...
                'Solution Mean is not within 1% Tolerance');
        end

        function TimeVaryingPoissonSpeedSuccess(testCase)
            % In this test, we check that the Time Varying 1D Poisson model
            % is solved in a reasonable amount of time.
            disp(['TV Poiss time = ',num2str(testCase.TvPoissSolution.time)]);
            testCase.verifyEqual(testCase.TvPoissSolution.time<0.3, true, ...
                'Time Varying Possion Solution Time is Slow');
        end        

        function ODEsolutionTV(testCase)
            % In this test, we check that the ODE Solution matches the
            % exact solution for the 1D Time Varying Poisson model and the
            % 2D Poisson Model.
            t = testCase.PoissTVODE.tSpan;
            tp = max(0,t-1);
            
            mn = testCase.PoissTVODE.parameters{1,2}/testCase.PoissTVODE.parameters{2,2}*...
                (1-exp(-testCase.PoissTVODE.parameters{2,2}*tp));
            
            Model = testCase.PoissTVODE;
            odeSoln = Model.solve;

            diff = abs(odeSoln.ode-mn');
            relDiff = sum(diff(mn>0)./mn(mn>0)')/length(odeSoln.ode);

            testCase.verifyEqual(relDiff<0.001, true, ...
                'ODE Solution is not within 1% Tolerance');
            
        end

    end
end