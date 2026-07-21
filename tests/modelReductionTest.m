classdef modelReductionTest < matlab.unittest.TestCase
    properties
        TwoDPoiss
        TwoDPoissSolution
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath(genpath('../src'));

            %% Test Case 3 - 2 Species Poisson Model
            tc.TwoDPoiss = SSIT;
            tc.TwoDPoiss.species = {'goff','gon','rna','prot'};
            tc.TwoDPoiss.initialCondition = [2;0;0;0];
            tc.TwoDPoiss.propensityFunctions = {'kon*goff*(1+sin(2*pi*t))';'koff*gon';'kr*gon';'gr*rna';'kp*rna';'gp*prot'};
            tc.TwoDPoiss.stoichiometry = [-1,1,0,0,0,0;1,-1,0,0,0,0;0,0,1,-1,0,0;0,0,0,0,1,-1];
            tc.TwoDPoiss.parameters = ({'kon',0.5';'koff',1;'kr',20;'gr',1;'kp',5;'gp',1});
            tc.TwoDPoiss.tSpan = linspace(0,5,201);
            tc.TwoDPoiss = tc.TwoDPoiss.formPropensitiesGeneral('TwoPoiss');

            [tc.TwoDPoissSolution,~,tc.TwoDPoiss] = tc.TwoDPoiss.solve(returnType='soln');
            tic
            [tc.TwoDPoissSolution,~,tc.TwoDPoiss] = tc.TwoDPoiss.solve(tc.TwoDPoissSolution.stateSpace,returnType='soln');
            tc.TwoDPoissSolution.time = toc;
     
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)

        % function NoTransform(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'No Transform';
        %     Model2.modelReductionOptions.reductionOrder = 25;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;
        % 
        %     Model2 = Model2.solve(fspSets.stateSpace);
        %     tic
        %     [fspSoln2,~,Model2] = Model2.solve(returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     NoTransform_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     NoTransformfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     NoTransformfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
        %     NoTransformfinalError = (NoTransformfinalError1+NoTransformfinalError2)/2;
        % 
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])
        %     % 
        %     % Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     % Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])
        % 
        %     tc.verifyEqual(NoTransformfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end

        % function LogLumpQSSA(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'Log Lump QSSA';
        %     Model2.modelReductionOptions.reductionOrder = 25;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;
        % 
        %     Model = Model2.solve(fspSets.stateSpace);
        %     tic
        %     [fspSoln2,~,Model2] = Model2.solve(returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     LogLumpQSSA_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     LogLumpQSSAfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     LogLumpQSSAfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
        %     LogLumpQSSAfinalError = (LogLumpQSSAfinalError1+LogLumpQSSAfinalError2)/2;
        % 
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2:5])
        %     % 
        %     % Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     % Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2:5])
        % 
        %     tc.verifyEqual(LogLumpQSSA_SpeedUp_Factor>1, true, ...
        %         'POD Reduction was not faster than original');
        %     tc.verifyEqual(LogLumpQSSAfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end

        % function EigenDecomposition(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'Eigen Decomposition';
        %     Model2.modelReductionOptions.reductionOrder = 25;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;
        % 
        %     Model = Model2.solve(fspSets.stateSpace);
        %     tic
        %     [fspSoln2,~,Model2] = Model2.solve(returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     EigenDecomposition_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     EigenDecompositionfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     EigenDecompositionfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
        %     EigenDecompositionfinalError = (EigenDecompositionfinalError1+EigenDecompositionfinalError2)/2;
        % 
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2:5])
        %     % 
        %     % Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     % Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2:5])
        % 
        %     tc.verifyEqual(EigenDecomposition_SpeedUp_Factor>1, true, ...
        %         'POD Reduction was not faster than original');
        %     tc.verifyEqual(EigenDecompositionfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end

        % function LinearStateLumping(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'Linear State Lumping';
        %     Model2.modelReductionOptions.reductionOrder = 30;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;
        % 
        %     Model2 = Model2.solve(fspSets.stateSpace);
        %     tic
        %     fspSoln2 = Model2.solve(returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     LinearStateLumping_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     LinearStateLumpingfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     LinearStateLumpingfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %     squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
        %     LinearStateLumpingfinalError = (LinearStateLumpingfinalError1+LinearStateLumpingfinalError2)/2;
        % 
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2:5])
        %     % 
        %     % Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     % Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2:5])
        % 
        %     tc.verifyEqual(LinearStateLumping_SpeedUp_Factor>1, true, ...
        %         'Reduction was not faster than original');
        %     tc.verifyEqual(LinearStateLumpingfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end

        % function LogarithmicStateLumping(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'Logarithmic State Lumping';
        %     Model2.modelReductionOptions.reductionOrder = 30;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;
        % 
        %     Model2 = Model2.solve(fspSets.stateSpace);
        %     tic
        %     fspSoln2 = Model2.solve(returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     LogarithmicStateLumping_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     LogarithmicStateLumpingfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     LogarithmicStateLumpingfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
        %     LogarithmicStateLumpingfinalError = (LogarithmicStateLumpingfinalError1+LogarithmicStateLumpingfinalError2)/2;
        % 
        %     tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2:5])
        % 
        %     Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2:5])
        % 
        %     tc.verifyEqual(LogarithmicStateLumping_SpeedUp_Factor>1, true, ...
        %         'Reduction was not faster than original');
        %     tc.verifyEqual(LogarithmicStateLumpingfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end

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
        %     fspSoln2 = Model2.solve(fspSets.stateSpace,returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     BalancedTruncation_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     BalancedTruncationfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %                  squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     BalancedTruncationfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %                  squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
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

            Model = Model2.solve(tc.TwoDPoissSolution.stateSpace);
            tic
            fspSoln2 = Model2.solve(tc.TwoDPoissSolution.stateSpace,returnType='soln');
            redModelSolveTime = toc;

            POD_SpeedUpFactor = tc.TwoDPoissSolution.time/redModelSolveTime
            PODfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
                squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
            PODfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
                squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
            PODfinalError = (PODfinalError1+PODfinalError2)/2;

            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2:5])

            Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2:5])

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

            Model2 = Model2.solve(tc.TwoDPoissSolution.stateSpace);
            tic
            fspSoln2 = Model2.solve(tc.TwoDPoissSolution.stateSpace,returnType='soln');
            redModelSolveTime = toc;

            POD_SpeedUpFactor = tc.TwoDPoissSolution.time/redModelSolveTime
            PODfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
                squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
            PODfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
            squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
            PODfinalError = (PODfinalError1+PODfinalError2)/2;

            % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
            % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])
            % 
            % Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
            % Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])

            tc.verifyEqual(POD_SpeedUpFactor>1, true, ...
                'POD Reduction was not faster than original');
            tc.verifyEqual(PODfinalError<0.01, true, ...
                'POD Reduction was not accurate enough');
        end

        % function QSSA(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'QSSA';
        %     Model2.modelReductionOptions.qssaSpecies = 2;
        %     Model2.modelReductionOptions.reductionOrder = 20;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices;
        % 
        %     Model2 = Model2.solve(fspSets.stateSpace);
        %     tic
        %     [fspSoln2,~,Model2] = Model2.solve(returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     QSSA_SpeedUpFactor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     QSSAfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     QSSAfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
        %     QSSAfinalError = (QSSAfinalError1+QSSAfinalError2)/2;
        % 
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])
        %     % 
        %     % Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     % Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])
        % 
        %     tc.verifyEqual(QSSA_SpeedUpFactor>1, true, ...
        %         'QSSA Reduction was not faster than original');
        %     tc.verifyEqual(QSSAfinalError<0.01, true, ...
        %         'QSSA Reduction was not accurate enough');
        % end

        % function DynamicModeDecomposition(tc)
        %     % In this test, we check that the 2D Poisson model generates a
        %     % solution with the correct means versus time.
        %     Model2 = tc.TwoDPoiss;
        %     Model2.modelReductionOptions.useModReduction = true;
        %     Model2.fspOptions.fspTol = inf;
        %     Model2.modelReductionOptions.reductionType = 'Dynamic Mode Decomposition';
        %     Model2.modelReductionOptions.reductionOrder = 25;
        %     [Model2,fspSets] = Model2.computeModelReductionTransformMatrices(tc.TwoDPoissSolution);
        % 
        %     Model2 = Model2.solve(fspSets.stateSpace);
        %     tic
        %     [fspSoln2,~,Model2] = Model2.solve(returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     DynamicModeDecomposition_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     DynamicModeDecompositionfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     DynamicModeDecompositionfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %         squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
        %     DynamicModeDecompositionfinalError = (DynamicModeDecompositionfinalError1+DynamicModeDecompositionfinalError2)/2;
        % 
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'meansAndDevs',[1:10:201],[],1)
        %     % tc.TwoDPoiss.makePlot(tc.TwoDPoissSolution,'marginals',[51,101,151,201],[],[2,3])
        %     % 
        %     % Model2.makePlot(fspSoln2,'meansAndDevs',[1:10:201],[],1)
        %     % Model2.makePlot(fspSoln2,'marginals',[51,101,152,201],[],[2,3])
        % 
        %     tc.verifyEqual(DynamicModeDecomposition_SpeedUp_Factor>1, true, ...
        %         'Reduction was not faster than original');
        %     tc.verifyEqual(DynamicModeDecompositionfinalError<0.05, true, ...
        %         'Transform Reduction was not accurate enough');
        % end

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
        %     fspSoln2 = Model2.solve(fspSets.stateSpace,returnType='soln');
        %     redModelSolveTime = toc;
        % 
        %     RadialBasisFunctions_SpeedUp_Factor = tc.TwoDPoissSolution.time/redModelSolveTime
        %     RadialBasisFunctionsfinalError1 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,4]))) - ...

        %               squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,4])))));
        %     RadialBasisFunctionsfinalError2 = max(abs(squeeze(cumsum(sum(double(fspSoln2.fsp{end}.p.data),[1,2,3]))) - ...
        %                  squeeze(cumsum(sum(double(tc.TwoDPoissSolution.fsp{end}.p.data),[1,2,3])))));
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