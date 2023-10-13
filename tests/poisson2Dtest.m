classdef poisson2Dtest < matlab.unittest.TestCase
    properties
        TwoDPoiss
        TwoDPoissSolution
        TwoDPoissODE
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(testCase1)
            addpath('../CommandLine')

            %% Test Case 3 - 2 Species Poisson Model
            testCase1.TwoDPoiss = SSIT;
            testCase1.TwoDPoiss.species = {'rna1','rna2'};
            testCase1.TwoDPoiss.initialCondition = [0;0];
            testCase1.TwoDPoiss.propensityFunctions = {'kr1';'gr1*rna1';'kr2';'gr2*rna2'};
            testCase1.TwoDPoiss.stoichiometry = [1,-1,0,0;0,0,1,-1];
            testCase1.TwoDPoiss.parameters = ({'kr1',10;'gr1',1;'kr2',5;'gr2',1});
            testCase1.TwoDPoiss.tSpan = linspace(0,2,21);
            testCase1.TwoDPoiss = testCase1.TwoDPoiss.formPropensitiesGeneral('TwoPoiss');

            [testCase1.TwoDPoissSolution,testCase1.TwoDPoiss.fspOptions.bounds] = testCase1.TwoDPoiss.solve;
            tic
            [testCase1.TwoDPoissSolution,testCase1.TwoDPoiss.fspOptions.bounds] = testCase1.TwoDPoiss.solve(testCase1.TwoDPoissSolution.stateSpace);
            testCase1.TwoDPoissSolution.time = toc;

            %% ODE model for Two Poisson process
            testCase1.TwoDPoissODE = testCase1.TwoDPoiss;
            testCase1.TwoDPoissODE.solutionScheme = 'ode';
            testCase1.TwoDPoissODE = testCase1.TwoDPoissODE.formPropensitiesGeneral('TwoDPoissODE');  

         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        
        function TwoDPoissMeansSuccess(testCase)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            t = testCase.TwoDPoiss.tSpan;
            mn1 = testCase.TwoDPoiss.parameters{1,2}/testCase.TwoDPoiss.parameters{2,2}*...
                (1-exp(-testCase.TwoDPoiss.parameters{2,2}*t));
            mn2 = testCase.TwoDPoiss.parameters{3,2}/testCase.TwoDPoiss.parameters{4,2}*...
                (1-exp(-testCase.TwoDPoiss.parameters{4,2}*t));
            
            fspSoln = testCase.TwoDPoissSolution.fsp;
            fspMn1 = NaN*mn1;
            fspMn2 = NaN*mn2;
            for i = 1:length(fspMn1)
                sz = size(testCase.TwoDPoissSolution.fsp{i}.p.data);                
                fspMn1(i) = [0:sz(1)-1]*double(testCase.TwoDPoissSolution.fsp{i}.p.sumOver(2).data);
                fspMn2(i) = [0:sz(2)-1]*double(testCase.TwoDPoissSolution.fsp{i}.p.sumOver(1).data);
            end

            diff1 = abs(fspMn1-mn1);
            diff2 = abs(fspMn2-mn2);

            relDiff1 = sum(diff1(diff1>0)./mn1(diff1>0))/length(fspSoln);
            relDiff2 = sum(diff2(diff2>0)./mn2(diff2>0))/length(fspSoln);

            testCase.verifyEqual(max(relDiff1,relDiff2)<0.01, true, ...
                'Two D Poiosson Solution Mean is not within 1% Tolerance');
        end

        function TwoDPoissonSpeedSuccess(testCase)
            % In this test, we check that the 2D Poisson model is solved in
            % a reasonable amount of time.
            disp(['2D Poiss time = ',num2str(testCase.TwoDPoissSolution.time)]);
            testCase.verifyEqual(testCase.TwoDPoissSolution.time<0.3, true, ...
                'TwoD Possion Solution Time is Slow');
        end

        function HybridModel(testCase)
            % In this test, we check that the code correctly solves the
            % hybrid 2D Poisson Model, where the first species is solved
            % using ODE model and the second species is solved using the
            % discrete stochastic analysis.

            HybridModel = testCase.TwoDPoiss;
            HybridModel.useHybrid = true;
            HybridModel.hybridOptions.upstreamODEs = {'rna1'};
            HybridModel = HybridModel.formPropensitiesGeneral('Hybrid');
            [hybSoln, HybridModel.fspOptions.bounds] = HybridModel.solve;

            t = HybridModel.tSpan;
            mn1 = HybridModel.parameters{1,2}/HybridModel.parameters{2,2}*...
                (1-exp(-HybridModel.parameters{2,2}*t));
            mn2 = HybridModel.parameters{3,2}/HybridModel.parameters{4,2}*...
                (1-exp(-HybridModel.parameters{4,2}*t));
            
            fspSoln = hybSoln.fsp;
            fspMn1 = NaN*mn1;
            fspMn2 = NaN*mn2;
            for i = 1:length(fspMn1)
                sz = size(fspSoln{i}.p.data,1);                
                fspMn2(i) = [0:sz-1]*double(fspSoln{i}.p.data);
                fspMn1(i) = fspSoln{i}.upstreamODEs;
            end

            diff1 = abs(fspMn1-mn1);
            diff2 = abs(fspMn2-mn2);

            relDiff1 = sum(diff1(diff1>0)./mn1(diff1>0))/length(fspSoln);
            relDiff2 = sum(diff2(diff2>0)./mn2(diff2>0))/length(fspSoln);

            testCase.verifyEqual(max(relDiff1,relDiff2)<0.01, true, ...
                'Hybrid Poisson Solution Mean is not within 1% Tolerance');
        end

        function Plotting(testCase)
            % In this test, we check that the code successfully generates
            % fogures for the solutions of the 2D Poisson Model. If
            % successful, four figures should have been generated.
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'meansAndDevs',[],[],1);
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'marginals',[3:3:21],[],[2,3]);
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'joints',[3:3:21],[],4);
        end

        function ODEsolution(testCase)
            % In this test, we check that the ODE Solution matches the
            % exact solution for the 1D Time Varying Poisson model and the
            % 2D Poisson Model.
            % Compare to 2D Poisson starting at SS.

            mn1 = testCase.TwoDPoissODE.parameters{1,2}/testCase.TwoDPoissODE.parameters{2,2};
            mn2 = testCase.TwoDPoissODE.parameters{3,2}/testCase.TwoDPoissODE.parameters{4,2};
            
            Model = testCase.TwoDPoissODE;
            Model.fspOptions.initApproxSS = true;
            odeSoln = Model.solve;

            diff1 = abs(odeSoln.ode(:,1)-mn1);
            diff2 = abs(odeSoln.ode(:,2)-mn2);

            relDiff1 = sum(diff1(diff1>0)./mn1)/length(diff1);
            relDiff2 = sum(diff2(diff2>0)./mn2)/length(diff2);

            testCase.verifyEqual(max([relDiff1,relDiff2])<0.001, true, ...
                'ODE Solution is not within 1% Tolerance');
            
        end
   
    end
end