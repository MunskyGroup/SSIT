classdef modelReductionTest < matlab.unittest.TestCase
    properties
        TwoDPoiss
        TwoDPoissSolution
        TwoDPoissODE
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath('../CommandLine')

            %% Test Case 3 - 2 Species Poisson Model
            tc.TwoDPoiss = SSIT;
            tc.TwoDPoiss.species = {'rna1','rna2'};
            tc.TwoDPoiss.initialCondition = [0;0];
            tc.TwoDPoiss.propensityFunctions = {'kr1';'gr1*rna1';'kr2';'gr2*rna2'};
            tc.TwoDPoiss.stoichiometry = [1,-1,0,0;0,0,1,-1];
            tc.TwoDPoiss.parameters = ({'kr1',100;'gr1',1;'kr2',80;'gr2',1});
            tc.TwoDPoiss.tSpan = linspace(0,5,201);
            tc.TwoDPoiss = tc.TwoDPoiss.formPropensitiesGeneral('TwoPoiss');

            [tc.TwoDPoissSolution,tc.TwoDPoiss.fspOptions.bounds] = tc.TwoDPoiss.solve;
            tic
            [tc.TwoDPoissSolution,tc.TwoDPoiss.fspOptions.bounds] = tc.TwoDPoiss.solve(tc.TwoDPoissSolution.stateSpace);
            tc.TwoDPoissSolution.time = toc;

            %% ODE model for Two Poisson process
            tc.TwoDPoissODE = tc.TwoDPoiss;
            tc.TwoDPoissODE.solutionScheme = 'ode';
            tc.TwoDPoissODE = tc.TwoDPoissODE.formPropensitiesGeneral('TwoDPoissODE');  

         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        
        function POD(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'Proper Orthogonal Decomposition';
            Model2.modelReductionOptions.reductionOrder = 25;
            Model2 = Model2.computeModelReductionTransformMatrices(tc.TwoDPoissSolution);

            tic
            fspSoln2 = Model2.solve(tc.TwoDPoissSolution.stateSpace);
            redModelSolveTime = toc;

            POD_SpeedUpFactor = tc.TwoDPoissSolution.time/redModelSolveTime
            PODfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            PODfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            PODfinalError = (PODfinalError1+PODfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(POD_SpeedUpFactor>1, true, ...
                 'POD Reduction was not faster than original');
            tc.verifyEqual(PODfinalError<0.01, true, ...
                 'POD Reduction was not accurate enough');
        end

        function QSSA(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'QSSA';
            Model2.modelReductionOptions.qssaSpecies = 1;
            Model2.modelReductionOptions.reductionOrder = 20;
            Model2 = Model2.computeModelReductionTransformMatrices(tc.TwoDPoissSolution);

            tic
            fspSoln2 = Model2.solve(tc.TwoDPoissSolution.stateSpace);
            redModelSolveTime = toc;

            QSSA_SpeedUpFactor = tc.TwoDPoissSolution.time/redModelSolveTime
            QSSAfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            QSSAfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            QSSAfinalError = (QSSAfinalError1+QSSAfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(QSSA_SpeedUpFactor>1, true, ...
                 'QSSA Reduction was not faster than original');
            tc.verifyEqual(QSSAfinalError<0.01, true, ...
                 'QSSA Reduction was not accurate enough');
        end

        function LogarithmicBinning(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'Log Lump QSSA';
            Model2.modelReductionOptions.reductionOrder = 25;
            Model2 = Model2.computeModelReductionTransformMatrices(tc.TwoDPoissSolution);

            tic
            fspSoln2 = Model2.solve(tc.TwoDPoissSolution.stateSpace);
            redModelSolveTime = toc;

            LogLump_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
            LogLumpfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            LogLumpfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            LogLumpfinalError = (LogLumpfinalError1+LogLumpfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(LogLump_SpeedUp_Factor>1, true, ...
                 'LogLump QSSA Reduction was not faster than original');
            tc.verifyEqual(LogLumpfinalError<0.05, true, ...
                 'LogLump QSSA Reduction was not accurate enough');
        end
    end
end