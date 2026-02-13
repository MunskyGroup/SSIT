classdef modelReductionTest < matlab.unittest.TestCase
    properties
        TwoDPoiss
        TwoDPoissSolution
        TwoDPoissODE
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath(genpath('../src'));

            %% Test Case 3 - 2 Species Poisson Model
            tc.TwoDPoiss = SSIT;
            tc.TwoDPoiss.species = {'rna1','rna2'};
            tc.TwoDPoiss.initialCondition = [0;0];
            tc.TwoDPoiss.propensityFunctions = {'kr1';'gr1*rna1';'kr2';'gr2*rna2'};
            tc.TwoDPoiss.stoichiometry = [1,-1,0,0;0,0,1,-1];
            tc.TwoDPoiss.parameters = ({'kr1',100;'gr1',3;'kr2',80;'gr2',3});
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

        function NoTransform(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'No Transform';
            Model2.modelReductionOptions.reductionOrder = 25;
            [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;

            tic
            fspSoln2 = Model2.solve(fspSets.stateSpace);
            redModelSolveTime = toc;

            NoTransform_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
            NoTransformfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            NoTransformfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            NoTransformfinalError = (NoTransformfinalError1+NoTransformfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(NoTransformfinalError<0.05, true, ...
                'Transform Reduction was not accurate enough');
        end

        function LogLumpQSSA(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'Log Lump QSSA';
            Model2.modelReductionOptions.reductionOrder = 25;
            [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;

            tic
            fspSoln2 = Model2.solve(fspSets.stateSpace);
            redModelSolveTime = toc;

            LogLumpQSSA_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
            LogLumpQSSAfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            LogLumpQSSAfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            LogLumpQSSAfinalError = (LogLumpQSSAfinalError1+LogLumpQSSAfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(LogLumpQSSA_SpeedUp_Factor>1, true, ...
                'POD Reduction was not faster than original');
            tc.verifyEqual(LogLumpQSSAfinalError<0.05, true, ...
                'Transform Reduction was not accurate enough');
        end

        function EigenDecomposition(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'Eigen Decomposition';
            Model2.modelReductionOptions.reductionOrder = 25;
            [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;

            tic
            fspSoln2 = Model2.solve(fspSets.stateSpace);
            redModelSolveTime = toc;

            EigenDecomposition_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
            EigenDecompositionfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            EigenDecompositionfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            EigenDecompositionfinalError = (EigenDecompositionfinalError1+EigenDecompositionfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(EigenDecomposition_SpeedUp_Factor>1, true, ...
                'POD Reduction was not faster than original');
            tc.verifyEqual(EigenDecompositionfinalError<0.05, true, ...
                'Transform Reduction was not accurate enough');
        end

        function LinearStateLumping(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'Linear State Lumping';
            Model2.modelReductionOptions.reductionOrder = 25;
            [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;

            tic
            fspSoln2 = Model2.solve(fspSets.stateSpace);
            redModelSolveTime = toc;

            LinearStateLumping_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
            LinearStateLumpingfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            LinearStateLumpingfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            LinearStateLumpingfinalError = (LinearStateLumpingfinalError1+LinearStateLumpingfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(LinearStateLumping_SpeedUp_Factor>1, true, ...
                'Reduction was not faster than original');
            tc.verifyEqual(LinearStateLumpingfinalError<0.05, true, ...
                'Transform Reduction was not accurate enough');
        end

        function LogarithmicStateLumping(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'Logarithmic State Lumping';
            Model2.modelReductionOptions.reductionOrder = 25;
            [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;

            tic
            fspSoln2 = Model2.solve(fspSets.stateSpace);
            redModelSolveTime = toc;

            LogarithmicStateLumping_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
            LogarithmicStateLumpingfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            LogarithmicStateLumpingfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            LogarithmicStateLumpingfinalError = (LogarithmicStateLumpingfinalError1+LogarithmicStateLumpingfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(LogarithmicStateLumping_SpeedUp_Factor>1, true, ...
                'Reduction was not faster than original');
            tc.verifyEqual(LogarithmicStateLumpingfinalError<0.05, true, ...
                'Transform Reduction was not accurate enough');
        end

        %% function BalancedModelTruncation(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'Balanced Model Truncation (HSV)';
        %     Model2.modelReductionOptions.reductionOrder = 25;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices(tc.TwoDPoissSolution);
        % 
        %     tic
        %     fspSoln2 = Model2.solve(fspSets.stateSpace);
        %     redModelSolveTime = toc;
        % 
        %     BalancedTruncation_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     BalancedTruncationfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
        %     BalancedTruncationfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
        %     BalancedTruncationfinalError = (BalancedTruncationfinalError1+BalancedTruncationfinalError2)/2;
        % 
        %     tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])
        % 
        %     Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])
        % 
        %     tc.verifyEqual(BalancedTruncation_SpeedUp_Factor>1, true, ...
        %         'Reduction was not faster than original');
        %     tc.verifyEqual(BalancedTruncationfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end
%%
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

        function POD2nd(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'POD 2nd';
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
            [Model2,fspSet] = Model2.computeModelReductionTransformMatrices;

            tic
            fspSoln2 = Model2.solve(fspSet.stateSpace);
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

        function DynamicModeDecomposition(tc)
            % In this test, we check that the 2D Poisson model generates a
            % solution with the correct means versus time.
            Model2 = tc.TwoDPoiss;
            Model2.modelReductionOptions.useModReduction = true;
            Model2.fspOptions.fspTol = inf;
            Model2.modelReductionOptions.reductionType = 'Dynamic Mode Decomposition';
            Model2.modelReductionOptions.reductionOrder = 25;
            [Model2,fspSets] = Model2.computeModelReductionTransformMatrices(tc.TwoDPoissSolution);

            tic
            fspSoln2 = Model2.solve(fspSets.stateSpace);
            redModelSolveTime = toc;

            DynamicModeDecomposition_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
            DynamicModeDecompositionfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
            DynamicModeDecompositionfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
            DynamicModeDecompositionfinalError = (DynamicModeDecompositionfinalError1+DynamicModeDecompositionfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(DynamicModeDecomposition_SpeedUp_Factor>1, true, ...
                'Reduction was not faster than original');
            tc.verifyEqual(DynamicModeDecompositionfinalError<0.05, true, ...
                'Transform Reduction was not accurate enough');
        end

        %% function RadialBasisFunctions(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'Radial Basis Functions';
        %     Model2.modelReductionOptions.reductionOrder = 200;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices(tc.TwoDPoissSolution);
        % 
        %     tic
        %     fspSoln2 = Model2.solve(fspSets.stateSpace);
        %     redModelSolveTime = toc;
        % 
        %     RadialBasisFunctions_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     RadialBasisFunctionsfinalError1 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),1))
        %     RadialBasisFunctionsfinalError2 = max(sum(abs((double(fspSoln2.fsp{end}.p.data - tc.TwoDPoissSolution.fsp{end}.p.data))),2))
        %     RadialBasisFunctionsfinalError = (RadialBasisFunctionsfinalError1+RadialBasisFunctionsfinalError2)/2;
        % 
        %     tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])
        % 
        %     Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])
        % 
        %     tc.verifyEqual(RadialBasisFunctions_SpeedUp_Factor>1, true, ...
        %         'Reduction was not faster than original');
        %     tc.verifyEqual(RadialBasisFunctionsfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end

    end
end