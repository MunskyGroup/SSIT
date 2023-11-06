classdef ssitTests < matlab.unittest.TestCase
    properties
        Poiss
        PoissSolution
        PoissODE
        TvPoiss
        TvPoissSolution
        PoissTVODE
        TwoDPoiss
        TwoDPoissSolution
        TwoDPoissODE
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(testCase1)
            addpath('../CommandLine')
            %% Test Case 1 - a simple Poisson model
            testCase1.Poiss = SSIT;
            testCase1.Poiss.species = {'rna'};
            testCase1.Poiss.initialCondition = 0;
            testCase1.Poiss.propensityFunctions = {'kr';'gr*rna'};
            testCase1.Poiss.stoichiometry = [1,-1];
            testCase1.Poiss.parameters = ({'kr',10;'gr',1});
            testCase1.Poiss.tSpan = linspace(0,2,21);
            testCase1.Poiss.fspOptions.fspTol = 1e-5;
            testCase1.Poiss = testCase1.Poiss.formPropensitiesGeneral('Poiss');

            [testCase1.PoissSolution,testCase1.Poiss.fspOptions.bounds] = testCase1.Poiss.solve;
            tic
            [testCase1.PoissSolution,testCase1.Poiss.fspOptions.bounds] = testCase1.Poiss.solve(testCase1.PoissSolution.stateSpace);
            testCase1.PoissSolution.time = toc;

            testCase1.Poiss = testCase1.Poiss.loadData('testData.csv',{'rna','exp1_s1'});


            %% ODE model for Poisson process
            testCase1.PoissODE = testCase1.Poiss;
            testCase1.PoissODE.solutionScheme = 'ode';
            testCase1.PoissODE = testCase1.PoissODE.formPropensitiesGeneral('PoissODE');  

            %% Test Case 2 - Poisson model with time-varying production
            testCase1.TvPoiss = testCase1.Poiss;
            testCase1.TvPoiss.propensitiesGeneral = []; % Processed propensity functions for use in solvers
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
        function ModelCreation(testCase)
            % In this trivial test, we check that the SSIT is set up with
            % the right names for the 'rna' species.
            nm = testCase.Poiss.species;
            testCase.verifyEqual(nm{1}, 'rna', ...
                'Species name is incorrect');
        end

        function FspConverged(testCase)
            % In this test, we check tha tthe FSP solution exits with an
            % appropriate FSP tolerance value.
            final = testCase.PoissSolution.fsp{end}.p.sum;
            tst = (1-final)<=testCase.Poiss.fspOptions.fspTol;
            testCase.verifyEqual(tst, true, ...
                'Final FSP is not within tolerance');
        end

        function PoissonMeans(testCase)
            % In this test, we check that the Time invariant 1D Poisson
            % model generates a solution with the correct means versus
            % time. 
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
            disp(['Poiss time = ',num2str(testCase.PoissSolution.time)]);
            testCase.verifyEqual(testCase.PoissSolution.time<0.2, true, ...
                'Possion Solution Time is Slow');
        end

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

        function EscapeTimeCalculation(testCase)
            % In this test, we check that the code correctly calcuates the
            % escape times until N births occur in the pure birth model.
            % Here, we compare to the exact solution:
            % P(t) = gamma(t|N,1/k)

            escapeModel = testCase.Poiss;
            escapeModel.parameters{2,2} = 0;
            t = escapeModel.tSpan;
            k = escapeModel.parameters{1,2};
            n = 10;
            escapeModel.fspOptions.escapeSinks.f = {'x1'};
            escapeModel.fspOptions.escapeSinks.b = n-0.1;
            % escapeModel = escapeModel.formPropensitiesGeneral;

            [escapeSoln,escapeModel.fspOptions.bounds] = escapeModel.solve;
            
            exact = cdf('gamma',t,n,1/k);

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
            % In this test, we check that the code successfully generates
            % fogures for the solutions of the 2D Poisson Model. If
            % successful, four figures should have been generated.
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'meansAndDevs',[],[],1);
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'marginals',[3:3:21],[],[2,3]);
            testCase.TwoDPoiss.makePlot(testCase.TwoDPoissSolution,'joints',[3:3:21],[],4);
        end

        function FspDataGeneration(testCase)
            % In this test, we check that the code generates and saves data
            % generated using the Poisson model.
            delete 'testData.csv'
            testCase.Poiss.ssaOptions.nSimsPerExpt = 1000;
            testCase.Poiss.ssaOptions.Nexp = 1;
            testCase.Poiss.sampleDataFromFSP(testCase.PoissSolution,'testData.csv')
            testCase.verifyEqual(exist('testData.csv','file'), 2, ...
                'FSP Data Not Generated');
            
        end

        function SsaDataGeneration(testCase)
            % In this test, we check that the code generates and saves data
            % generated using the Poisson model.
            delete 'testDataSSA.csv'
            testCase.Poiss.ssaOptions.nSimsPerExpt = 1000;
            testCase.Poiss.ssaOptions.Nexp = 1;
            testCase.Poiss.solutionScheme = 'SSA';
            testCase.Poiss.solve(testCase.PoissSolution,'testDataSSA.csv')
            testCase.verifyEqual(exist('testDataSSA.csv','file'), 2, ...
                'SSA Data Not Generated');            
        end

        function DataLoading(testCase)
            % In this test, we check that the code correctly loads the data
            % geerated using the Poisson model. 
            % Note - this will only work if the code successfully generated
            % the data in the first place.
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
            % In this test we compare the computed log-likelihood to the
            % exact solution for the Poisson Model:
            % lam(t) = k/g*(1-exp(-g*t));
            % logL = prod_n [Poisson(n|lam(t))]
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

        function ComputingSensitivities(testCase)
            % In this test, we compare the poisson model sensitivites to
            % the parameters 'k' and 'g' using the exact solution:
            % lam(t) = k/g*(1-exp(-g*t));
            % Sk(t) = dP(n|t)/dk = 1/g*(1-exp(-g*t)) * Poiss(n|lam(t)) * (n/lam - 1)
            % Sg(t) = dP(n|t)/dg = (-k/g^2*(1-exp(-g*t))+k*t/g*exp(-g*t)) * Poiss(n|lam(t)) * (n/lam - 1)

            Model = testCase.Poiss;
            % Model = Model.formPropensitiesGeneral;
            Model.solutionScheme = 'FSP';
            [~,Model.fspOptions.bounds] = Model.solve;

            
            Model.solutionScheme = 'fspSens';
            Model.fspOptions.fspTol = 1e-6;            

            SensSoln = Model.solve;
            
            t = Model.tSpan;
            k = Model.parameters{1,2};
            g = Model.parameters{2,2};
            lam = k/g*(1-exp(-g*t))';

            fspSensSoln = SensSoln.sens.data;
            for i = 2:length(fspSensSoln)
                n = size(fspSensSoln{i}.S(1).data);
                analytical1 = pdf('poiss',(0:n-1),lam(i))'.*([0:n-1]'/lam(i) - 1).*...
                    (1/g)*(1-exp(-g*t(i)));
                fspSens1 = double(fspSensSoln{i}.S(1).data);
                diff1(i) = sum(abs(analytical1-fspSens1));
                
                analytical2 = pdf('poiss',(0:n-1),lam(i))'.*([0:n-1]'/lam(i) - 1).*...
                    ((-k/g^2)*(1-exp(-g*t(i)))+k/g*t(i)*exp(-g*t(i)));
                fspSens2 = double(SensSoln.sens.data{i}.S(2).data);
                diff2(i) = sum(abs(analytical2-fspSens2));
                
            end

            Model.makePlot(SensSoln,'marginals',21,[],[5]);

            testCase.verifyEqual(max(diff1+diff2)<0.001, true, ...
                'Sensitivity Calculation is not within 0.1% Tolerance');            

        end

        function TestFIM(testCase)
            % In this test case, we check that the FIM calculation using
            % the FSP matches to the analytical expresion for the Poisson
            % model:
            % lam(t) = k/g*(1-exp(-g*t));
            % FIM_kk = (dlam/dk)^2 * FIM_lamlam = lam/k^2
            % FIM_gg = (dlam/dg)^2 * FIM_lamlam = (-lam/g + k*t/g*exp(-g*t))
            % FIM_kg = (dlam/dg)(dlam/dk) * FIM_lamlam = 
            %                                 (-lam/g + k*t/g*exp(-g*t))/k
            Model = testCase.Poiss;
            Model.solutionScheme = 'fspSens';
            Model.fspOptions.fspTol = 1e-6;
            SensSoln = Model.solve;
            fspFIM = Model.computeFIM(SensSoln.sens);
            t = Model.tSpan;
            k = Model.parameters{1,2};
            g = Model.parameters{2,2};
            lam = k/g*(1-exp(-g*t))';
            for i = 2:length(t)
                exactIk(i) = lam(i)/k^2;
                fspIk(i) = fspFIM{i}(1,1);
                exactIg(i) = (-lam(i)/g + k*t(i)/g*exp(-g*t(i)))^2/lam(i);
                fspIg(i) = fspFIM{i}(2,2);
                exactIkg(i) = (-lam(i)/g + k*t(i)/g*exp(-g*t(i)))/k;
                fspIkg(i) = fspFIM{i}(1,2);
            end
                
            diff = max(abs(exactIk-fspIk)+abs(exactIg-fspIg)+abs(exactIkg-fspIkg));
            testCase.verifyEqual(diff<0.001, true, ...
                'FIM Calculation is not within 1e-4% Tolerance');            


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

            testCase.verifyEqual(max([relDiff,relDiff1,relDiff2])<0.001, true, ...
                'ODE Solution is not within 1% Tolerance');
            
        end
   
        function ComputeODELikelihood(testCase)
            % In this test we compare the computed ODE approximation for
            % the log-likelihood to the exact solution for the Poisson Model:
            % lam(t) = k/g*(1-exp(-g*t));
            % logL = prod_n [Poisson(n|lam(t))]
            Model = testCase.PoissODE;
            % Model.solutionScheme = 'ode';
            % Model = Model.formPropensitiesGeneral;            
            odeLogL = Model.computeLikelihoodODE;
            
            t = testCase.Poiss.tSpan;
            mn = testCase.Poiss.parameters{1,2}/testCase.Poiss.parameters{2,2}*...
                (1-exp(-testCase.Poiss.parameters{2,2}*t));

            logLExact = -1/2*1000*sum((testCase.Poiss.dataSet.mean-mn').^2);

            relDiff = abs((logLExact-odeLogL)/logLExact);

            testCase.verifyEqual(relDiff<0.001, true, ...
                'ODE Likelihood Calculation is not within 0.01% Tolerance');            
            
        end

        function FitODEModel(testCase)
            Model = testCase.PoissODE;
            % Model.solutionScheme = 'ode';
            % Model = Model.formPropensitiesGeneral;
            Model.parameters(:,2) = {10*rand;rand};
            fitOptions = optimset('Display','none','MaxIter',1000);
            fitOptions.SIG = [];
            Model.fittingOptions.modelVarsToFit = [1,2];
            for i=1:3
                fitPars = Model.maximizeLikelihood([],fitOptions);
                Model.parameters(:,2) = num2cell(fitPars);
            end

            relDiff = max(abs([testCase.Poiss.parameters{:,2}]-...
                [Model.parameters{:,2}])./[testCase.Poiss.parameters{:,2}]);

            testCase.verifyEqual(relDiff<0.05, true, ...
                'ODE Fit of Poisson Model is not within 5% of true values');            
        end 

        function FitUsingFSP(testCase)
            Model = testCase.Poiss;
            Model.solutionScheme = 'FSP';
            % Model = Model.formPropensitiesGeneral;
            Model.parameters(:,2) = {10*rand;rand};
            fitOptions = optimset('Display','none','MaxIter',1000);
            fitOptions.SIG = [];
            Model.fittingOptions.modelVarsToFit = [1,2];
            for i=1:3
                fitPars = Model.maximizeLikelihood([],fitOptions);
                Model.parameters(:,2) = num2cell(fitPars);
            end

            relDiff = max(abs([testCase.Poiss.parameters{:,2}]-...
                [Model.parameters{:,2}])./[testCase.Poiss.parameters{:,2}]);

            testCase.verifyEqual(relDiff<0.05, true, ...
                'ODE Fit of Poisson Model is not within 5% of true values');            
        end  
        
        function MetHastAndSampledFIM(testCase)
            Model = testCase.Poiss;
            Model.solutionScheme = 'FSP';
            % Model = Model.formPropensitiesGeneral;

            % run Metropolis Hastings
            MHFitOptions.thin=1;
            MHFitOptions.numberOfSamples=1000;
            MHFitOptions.burnIn=0;
            MHFitOptions.progress=true;
            MHFitOptions.useFIMforMetHast =true;
            MHFitOptions.CovFIMscale = 1.0;
            MHFitOptions.numChains = 1;
            MHFitOptions.saveFile = 'TMPMHChain.mat';
            Model.fittingOptions.modelVarsToFit = [1,2];
            [newPars,~,MHResults] = Model.maximizeLikelihood(...
                [], MHFitOptions, 'MetropolisHastings');
            Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(newPars);
            delete('TMPMHChain.mat')            

            % Compute FIM for subsampling of MH results.
            J = floor(linspace(500,1000,10));
            MHSamplesForFIM = exp(MHResults.mhSamples(J,:));
            fimResults = Model.computeFIM([],'lin',MHSamplesForFIM);

            % Find and plot total FIM for each parameter sample
            Nc = Model.dataSet.nCells;
            figure
            FIM = Model.totalFim(fimResults,Nc);
            Model.plotMHResults(MHResults,FIM)

            % Find optimal experiment design given parameters sets
            NcOptExperiment = Model.optimizeCellCounts(fimResults,sum(Nc),'Determinant');
            FIMOptExpt = Model.totalFim(fimResults,NcOptExperiment);
            Model.plotMHResults(MHResults,[FIM,FIMOptExpt])
            
        end  
    end
end