classdef ssitTests < matlab.unittest.TestCase
    properties
        Poiss
        PoissSolution
        TvPoiss
        TvPoissSolution
        TwoDPoiss
        TwoDPoissSolution
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(testCase1)
            addpath('CommandLine')
            %% Test Case 1 - a simple Poisson model
            testCase1.Poiss = SSIT;
            testCase1.Poiss.species = {'rna'};
            testCase1.Poiss.initialCondition = 0;
            testCase1.Poiss.propensityFunctions = {'kr';'gr*rna'};
            testCase1.Poiss.stoichiometry = [1,-1];
            testCase1.Poiss.parameters = ({'kr',10;'gr',1});
            testCase1.Poiss.tSpan = linspace(0,2,21);

            tic
            [testCase1.PoissSolution,testCase1.Poiss.fspOptions.bounds] = testCase1.Poiss.solve;
            testCase1.PoissSolution.time = toc;

            testCase1.Poiss = testCase1.Poiss.loadData('testData.csv',{'rna','exp1_s1'});

        end

        function createTestModel2(testCase2)
            %% Test Case 2 - Poisson model with time-varying production
            testCase2.TvPoiss = testCase2.Poiss;
            testCase2.TvPoiss.propensityFunctions = {'kr*Ig';'gr*rna'};
            testCase2.TvPoiss.inputExpressions = {'Ig','t>1'};
            
            tic
            [testCase2.TvPoissSolution,testCase2.TvPoiss.fspOptions.bounds] = testCase2.TvPoiss.solve;
            testCase2.TvPoissSolution.time = toc;
        end

        function createTestModel3(testCase3)
            %% Test Case 3 - 2 Species Poisson Model
            testCase3.TwoDPoiss = SSIT;
            testCase3.TwoDPoiss.species = {'rna1','rna2'};
            testCase3.TwoDPoiss.initialCondition = [0;0];
            testCase3.TwoDPoiss.propensityFunctions = {'kr1';'gr1*rna1';'kr2';'gr2*rna2'};
            testCase3.TwoDPoiss.stoichiometry = [1,-1,0,0;0,0,1,-1];
            testCase3.TwoDPoiss.parameters = ({'kr1',10;'gr1',1;'kr2',5;'gr2',1});
            testCase3.TwoDPoiss.tSpan = linspace(0,2,21);
            
            tic
            [testCase3.TwoDPoissSolution,testCase3.TwoDPoiss.fspOptions.bounds] = testCase3.TwoDPoiss.solve;
            testCase3.TwoDPoissSolution.time = toc;
            
         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function ModelCreation(testCase)
            nm = testCase.Poiss.species;
            testCase.verifyEqual(nm{1}, 'rna', ...
                'Species name is incorrect');
        end

        function FspConverged(testCase)
            final = testCase.PoissSolution.fsp{end}.p.sum;
            tst = (1-final)<=testCase.Poiss.fspOptions.fspTol;
            testCase.verifyEqual(tst, true, ...
                'Final FSP is not within tolerance');
        end

        function PoissonMeans(testCase)
            t = testCase.Poiss.tSpan;
            mn = testCase.Poiss.parameters{1,2}/testCase.Poiss.parameters{2,2}*...
                (1-exp(-testCase.Poiss.parameters{2,2}*t));
            
            fspSoln = testCase.PoissSolution.fsp;
            fspMn = NaN*mn;
            for i = 1:length(fspSoln)
                n = size(testCase.PoissSolution.fsp{i}.p.data,1);
                fspMn(i) = [0:n-1]*double(testCase.PoissSolution.fsp{i}.p.data);
            end

            diff = abs(fspMn-mn);
            relDiff = sum(diff(diff>0)./mn(diff>0))/length(fspSoln);

            testCase.verifyEqual(relDiff<0.01, true, ...
                'Solution Mean is not within 1% Tolerance');
        end

        function PoissonSpeed(testCase)
            testCase.verifyEqual(testCase.PoissSolution.time<0.2, true, ...
                'Possion Solution Time is Slow');
        end

        function TimeVaryingPoissonMeansSuccess(testCase)
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
            testCase.verifyEqual(testCase.TvPoissSolution.time<0.4, true, ...
                'Time Varying Possion Solution Time is Slow');
        end
        
        function TwoDPoissMeansSuccess(testCase)
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
            testCase.verifyEqual(testCase.TwoDPoissSolution.time<0.3, true, ...
                'TwoD Possion Solution Time is Slow');
        end

        function HybridModel(testCase)

            HybridModel = testCase.TwoDPoiss;
            HybridModel.useHybrid = true;
            HybridModel.hybridOptions.upstreamODEs = {'rna1'};
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

        function EscapeTimeCalculation(testCase)

            escapeModel = testCase.Poiss;
            escapeModel.parameters{2,2} = 0;
            escapeModel.fspOptions.escapeSinks.f = {'x1'};
            escapeModel.fspOptions.escapeSinks.b = 10;
            [escapeSoln,escapeModel.fspOptions.bounds] = escapeModel.solve;
            
            t = escapeModel.tSpan;
            a = escapeModel.parameters{1,2};
            n = 10;

            exact = cdf('gamma',t,n+1,1/a);

            fspSoln = escapeSoln.fsp;
            escapeProbFSP = NaN*exact;
            for i = 1:length(escapeProbFSP)
                escapeProbFSP(i) = fspSoln{i}.escapeProbs;
            end
          
            diff1 = abs(escapeProbFSP-exact);

            relDiff1 = sum(diff1(diff1>0)./exact(diff1>0))/length(escapeProbFSP);

            testCase.verifyEqual(relDiff1<0.01, true, ...
                'Pure Birth Escape Time CDF is not within 1% Tolerance');
           
        end

        function Plotting(testCase)
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'meansAndDevs',[],[],1);
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'marginals',[3:3:21],[],[2,3]);
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'joints',[3:3:21],[],[4]);
        end

        function FspDataGeneration(testCase)
            delete 'testData.csv'
            testCase.Poiss.ssaOptions.nSimsPerExpt = 1000;
            testCase.Poiss.ssaOptions.Nexp = 1;
            testCase.Poiss.sampleDataFromFSP(testCase.PoissSolution,'testData.csv')
            testCase.verifyEqual(exist('testData.csv','file'), 2, ...
                'FSP Data Not Generated');
            
        end

        function DataLoading(testCase)
            Z = double(testCase.Poiss.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor);
            n = sum(Z(1,:));
            ssaMn = Z*[0:size(Z,2)-1]'/n;
                        
            t = testCase.Poiss.tSpan;
            mn = testCase.Poiss.parameters{1,2}/testCase.Poiss.parameters{2,2}*...
                (1-exp(-testCase.Poiss.parameters{2,2}*t));
 

            diff = abs(ssaMn-mn');
            relDiff = sum(diff(diff>0)./mn(diff>0))/length(ssaMn);

            testCase.verifyEqual(max(relDiff)<0.1, true, ...
                'Generated Data may be inaccurate ');
            
        end

        function ComputingLikelihood(testCase)
            fspLogL = testCase.Poiss.computeLikelihood;
            
            t = [testCase.Poiss.dataSet.DATA{:,1}];
            x = [testCase.Poiss.dataSet.DATA{:,2}];
            mn = testCase.Poiss.parameters{1,2}/testCase.Poiss.parameters{2,2}*...
                (1-exp(-testCase.Poiss.parameters{2,2}*t));

            logLExact = sum(log(pdf('poiss',x,mn)));

            relDiff = abs((logLExact-fspLogL)/logLExact);

            testCase.verifyEqual(relDiff<0.0001, true, ...
                'Likelihood Calculation is not within 0.01% Tolerance');            
        end

    end
end