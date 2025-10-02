classdef poissonTest < matlab.unittest.TestCase
    properties
        Poiss
        PoissSolution
        PoissODE
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(testCase1)
            addpath(genpath('../src'));
            %% Test Case 1 - a simple Poisson model
            testCase1.Poiss = SSIT;
            testCase1.Poiss.species = {'rna'};
            testCase1.Poiss.initialCondition = 0;
            testCase1.Poiss.propensityFunctions = {'kr';'gr*rna'};
            testCase1.Poiss.stoichiometry = [1,-1];
            testCase1.Poiss.parameters = ({'kr',10;'gr',1});
            testCase1.Poiss.tSpan = linspace(0,2,21);
            testCase1.Poiss.fspOptions.fspTol = 1e-5;
            testCase1.Poiss = testCase1.Poiss.formPropensitiesGeneral('Poiss',true);

            [testCase1.PoissSolution, testCase1.Poiss.fspOptions.bounds, testCase1.Poiss] = testCase1.Poiss.solve;
            tic
            [testCase1.PoissSolution,testCase1.Poiss.fspOptions.bounds] = testCase1.Poiss.solve(testCase1.PoissSolution.stateSpace);
            testCase1.PoissSolution.time = toc;

            delete 'testData.csv'
            testCase1.Poiss.ssaOptions.nSimsPerExpt = 1000;
            testCase1.Poiss.ssaOptions.Nexp = 1;
            testCase1.Poiss.sampleDataFromFSP(testCase1.PoissSolution,'testData.csv');

            testCase1.Poiss = testCase1.Poiss.loadData('testData.csv',{'rna','exp1_s1'});

            %% ODE model for Poisson process
            testCase1.PoissODE = testCase1.Poiss;
            testCase1.PoissODE.solutionScheme = 'ode';
            testCase1.PoissODE = testCase1.PoissODE.formPropensitiesGeneral('PoissODE');  

         end  
    end

    methods (TestMethodSetup)
        
    end

    methods (Test)
        function ModelCreation(testCase)
            % In this trivial test, we check that the SSIT is set up with
            % the right names for the 'rna' species.
            nm = testCase.Poiss.species;
            testCase.verifyEqual(nm{1}, 'rna', ...
                'Species name is incorrect');
        end

        function AddReaction(testCase)
            % In this trivial test, we check that the SSIT is set up with
            % the right names for the 'rna' species.
            F = SSIT('Empty');
            newRxn(1).propensity = 'kr + kr1*x1';
            newRxn(1).stoichiometry = {'x1',1};
            newRxn(1).parameters = {'kr',2;'kr1',0.01};
            newRxn(2).propensity = 'g*x1';
            newRxn(2).stoichiometry = {'x1',-1};
            newRxn(2).parameters = {'g',0.1};
            F = F.addReaction(newRxn);

            testCase.verifyEqual(F.species{1}, 'x1', ...
                'FSP add reaction has error');
            testCase.verifyEqual(F.propensityFunctions{1}, 'kr + kr1*x1', ...
                'FSP add reaction has error');
            testCase.verifyEqual(F.stoichiometry, [1    -1], ...
                'FSP add reaction has error');
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

        function PoissonMoments(testCase)            
            % This test verifies that the moment closure approach is
            % working for the Poisson case.
            model = testCase.Poiss;

            t = testCase.Poiss.tSpan;
            mn = testCase.Poiss.parameters{1,2}/testCase.Poiss.parameters{2,2}*...
                (1-exp(-testCase.Poiss.parameters{2,2}*t));

            model.solutionScheme = 'ode';
            [~,~,model] = model.solve;

            model.solutionScheme = 'moments';
            [~,~,model] = model.solve;

            errMean = max(abs((model.Solutions.moments(1,:)-mn)./mn));

            var = model.Solutions.moments(2,:)-model.Solutions.moments(1,:).^2;
            errVar = max(abs((var-mn)./mn));

            testCase.verifyEqual(errMean+errVar<0.01, true, ...
                'Solution Mean and Variance not within 1% Tolerance');
            
        end

        function PoissonSpeed(testCase)
            disp(['Poiss time = ',num2str(testCase.PoissSolution.time)]);
            testCase.verifyEqual(testCase.PoissSolution.time<0.2, true, ...
                'Possion Solution Time is Slow');
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

        function FspDataGeneration(testCase)
            % In this test, we check that the code generates and saves data
            % generated using the Poisson model.
            delete 'testData.csv'
            testCase.Poiss.ssaOptions.nSimsPerExpt = 1000;
            testCase.Poiss.ssaOptions.Nexp = 1;
            testCase.Poiss.sampleDataFromFSP(testCase.PoissSolution,'testData.csv');
            testCase.verifyEqual(exist('testData.csv','file'), 2, ...
                'FSP Data Not Generated');
            
        end

        function SsaDataGeneration(testCase)
            % In this test, we check that the code generates and saves data
            % generated using the Poisson model.
            delete 'testDataSSA.csv'
            testCase.Poiss.ssaOptions.useParalel = true;
            testCase.Poiss.ssaOptions.nSimsPerExpt = 1000;
            testCase.Poiss.ssaOptions.Nexp = 1;
            testCase.Poiss.solutionScheme = 'SSA';
            testCase.Poiss.solve(testCase.PoissSolution,'testDataSSA.csv')
            testCase.verifyEqual(exist('testDataSSA.csv','file'), 2, ...
                'SSA Data Not Generated');            
        end

        function SSALossFunction(testCase)
            % In this test, we check that the code generates and saves data
            % generated using the Poisson model.
            testCase.Poiss.ssaOptions.useParalel = true;
            testCase.Poiss.ssaOptions.nSimsPerExpt = 1000;
            testCase.Poiss.ssaOptions.Nexp = 1;
            testCase.Poiss.solutionScheme = 'SSA';

            % Generate new data
            [~,~,testCase.Poiss] = testCase.Poiss.solve(testCase.PoissSolution,'testDataSSANEW.csv');
            
            % load new data and check that the loss function is now zero.
            testCase.Poiss = testCase.Poiss.loadData('testDataSSANEW.csv',{'rna','exp1_s1'});
            
            [lossFunction] = testCase.Poiss.computeLossFunctionSSA('cdf_one_norm',[],true,true);

            testCase.verifyEqual(lossFunction, 0, ...
                'SSA Loss Function Incorrectly Calculated');   

            % Test custom loss function
            customLossMeans = @(Hmod,Hdata)([0:length(Hmod)-1]*Hmod-[0:length(Hdata)-1]*Hdata).^2;
            [lossFunction] = testCase.Poiss.computeLossFunctionSSA(customLossMeans,[],true,true);

            testCase.verifyEqual(lossFunction, 0, ...
                'SSA Loss Function Incorrectly Calculated');

            % Test ABC for a short run.
            fitOptions.numberOfSamples=100;
            [~,~,Results] = testCase.Poiss.runABCsearch([],[],[],fitOptions);
            

        end

        function SsaDataGenerationGPU(testCase)
            % In this test, we check that the code generates and saves data
            % generated using the Poisson model.
            testCase.Poiss.ssaOptions.nSimsPerExpt = 100000;
            testCase.Poiss.ssaOptions.Nexp = 1;
            testCase.Poiss.solutionScheme = 'SSA';
            testCase.Poiss.tSpan = [0:10];

            try
                parpool
            catch
            end

            tic
            testCase.Poiss.ssaOptions.useParalel = false;
            testCase.Poiss.ssaOptions.useGPU = true;
            SSAGPU = testCase.Poiss.solve(testCase.PoissSolution);
            timeGPU = toc;
           
            tic
            testCase.Poiss.ssaOptions.useParallel = true;
            testCase.Poiss.ssaOptions.useGPU = false;
            SSACPUp = testCase.Poiss.solve(testCase.PoissSolution);
            timeCPUp = toc;

            tic
            testCase.Poiss.ssaOptions.useParallel = false;
            testCase.Poiss.ssaOptions.useGPU = false;
            SSACPUs = testCase.Poiss.solve(testCase.PoissSolution);
            timeCPUs = toc;    

            p = gcp("nocreate");
            Methods = {'1CPU+1GPU';append(num2str(p.NumWorkers),' CPUs');'1 CPU'};
            Times = [timeGPU;timeCPUp;timeCPUs];
            Means = [mean(SSAGPU.trajs(1,end,:));mean(SSACPUp.trajs(1,end,:));mean(SSACPUs.trajs(1,end,:))];
            Vars = [var(SSAGPU.trajs(1,end,:));var(SSACPUp.trajs(1,end,:));var(SSACPUs.trajs(1,end,:))];
            
            disp('GPU/CPU Performance for SSA - 1.1Mk sims')
            table(Methods,Times,Means,Vars)

            testCase.Poiss.ssaOptions.nSimsPerExpt = 1000000;
            tic
            testCase.Poiss.ssaOptions.useParalel = false;
            testCase.Poiss.ssaOptions.useGPU = true;
            SSAGPU = testCase.Poiss.solve(testCase.PoissSolution);
            timeGPU = toc;
           
            tic
            testCase.Poiss.ssaOptions.useParallel = true;
            testCase.Poiss.ssaOptions.useGPU = false;
            SSACPUp = testCase.Poiss.solve(testCase.PoissSolution);
            timeCPUp = toc;

            Methods = {'1CPU+1GPU';append(num2str(p.NumWorkers),' CPUs')};
            Times = [timeGPU;timeCPUp];
            Means = [mean(SSAGPU.trajs(1,end,:));mean(SSACPUp.trajs(1,end,:))];
            Vars = [var(SSAGPU.trajs(1,end,:));var(SSACPUp.trajs(1,end,:))];
            
            disp('GPU/CPU Performance for SSA - 11M sims')
            table(Methods,Times,Means,Vars)

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
            
            % Add additional times to FSP solution to test that it
            % correctly can filter thee out when computing the likelihood
            % values.
            testCase.Poiss.tSpan = [0:0.05:max(testCase.Poiss.tSpan)];
            
            % Update solution.
            [~,~,testCase.Poiss] = testCase.Poiss.solve;

            % Call to compute likelihood
            fspLogL = testCase.Poiss.computeLikelihood([],[],false,true);
            
            t = [testCase.Poiss.dataSet.DATA{:,1}];
            x = [testCase.Poiss.dataSet.DATA{:,2}];
            mn = testCase.Poiss.parameters{1,2}/testCase.Poiss.parameters{2,2}*...
                (1-exp(-testCase.Poiss.parameters{2,2}*t));

            logLExact = sum(log(pdf('poiss',x,mn)));

            relDiff = abs((logLExact-fspLogL)/logLExact);

            testCase.verifyEqual(relDiff<0.0001, true, ...
                'Likelihood Calculation is not within 0.01% Tolerance');            
        end

        function LikelihoodGradient(testCase)
            % This tests to make sure that the calculation for the gradient
            % of the loglikelihood function completes and is within 0.1% of
            % the solution found using the finite difference method.
            testCase.Poiss.solutionScheme = 'fspSens';
            [~,~,testCase.Poiss] = testCase.Poiss.solve;

            [fspLogL,gradient] = testCase.Poiss.computeLikelihood([],[],true,true);

            testCase.Poiss.solutionScheme = 'FSP';
            numPars = size(testCase.Poiss.parameters,1);
            gradLogLFiniteDiff = zeros(numPars,1);
            for i = 1:numPars
                TMPmodel = testCase.Poiss;
                delt = abs(TMPmodel.parameters{i,2})/1e6;
                TMPmodel.parameters{i,2} = TMPmodel.parameters{i,2} + delt;
                fspLogLPrime = TMPmodel.computeLikelihood;
                gradLogLFiniteDiff(i) = (fspLogLPrime-fspLogL)/delt;
            end

            maxGradientError = max(abs(gradient-gradLogLFiniteDiff)./gradLogLFiniteDiff);

            testCase.verifyEqual(maxGradientError<0.001, true, ...
                'Likelihood Gradient Calculation is not within 0.1% Tolerance');

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

            Model.sensOptions.solutionMethod = 'forward';
            SensSoln = Model.solve;
            Model.sensOptions.solutionMethod = 'finiteDifference';
            SensSoln2 = Model.solve;
            
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
            subplot(2,1,1)
            hold on
            plot(analytical1);
            subplot(2,1,2)
            hold on
            plot(analytical2);
            Model.makePlot(SensSoln2,'marginals',21,[],[5]);

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

        function TestFIMwPrior(testCase)
            % In this test case, we check that the FIM calculation using
            % the FSP matches to the analytical expresion for the Poisson
            % model when there is a appropriate Prior:
            Model = testCase.Poiss;
            Model.solutionScheme = 'fspSens';
            Model.fspOptions.fspTol = 1e-6;
            SensSoln = Model.solve;
            fspFIM = Model.computeFIM(SensSoln.sens);
            t = Model.tSpan;
            SIGprior = diag(rand(1,2));
            FIMwPrior = Model.evaluateExperiment(fspFIM, [1:length(t)],...
                SIGprior);
            FIMwoPrior = Model.evaluateExperiment(fspFIM,[1:length(t)]);
            diff = max(abs(FIMwPrior{1}-FIMwoPrior{1}-inv(SIGprior)),[],"all");
            testCase.verifyEqual(diff<1e-5, true, ...
                'FIM Prior Calculation is not within 1e-4% Tolerance');            

        end

        function ComputeODELikelihood(testCase)
            % In this test we compare the computed ODE approximation for
            % the log-likelihood to the exact solution for the Poisson Model:
            % lam(t) = k/g*(1-exp(-g*t));
            % logL = prod_n [Poisson(n|lam(t))]
            Model = testCase.PoissODE;
            
            % Skip first time point for fitting.
            Model.fittingOptions.timesToFit = Model.tSpan>0;
            odeLogL = Model.computeLikelihoodODE;
            
            t = Model.tSpan(2:end);
            mn = Model.parameters{1,2}/Model.parameters{2,2}*...
                (1-exp(-Model.parameters{2,2}*t));

            logLExact = -1/2*sum(Model.dataSet.nCells(2:end).*(Model.dataSet.mean(2:end)-mn').^2./Model.dataSet.var(2:end));

            relDiff = abs((logLExact-odeLogL)/logLExact);

            testCase.verifyEqual(relDiff<0.001, true, ...
                'ODE Likelihood Calculation is not within 0.01% Tolerance');            
            
        end

        function FitODEModel(testCase)
            Model = testCase.PoissODE;
            % Model.solutionScheme = 'ode';
            % Model = Model.formPropensitiesGeneral;
            Model.parameters(:,2) = {10*rand;rand};
            Model.fittingOptions.timesToFit = Model.tSpan>0;
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

            Model.makeFitPlot()
        end  
        
        function MetHastAndSampledFIM(testCase)
            Model = testCase.Poiss;
            Model.solutionScheme = 'FSP';

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
            NcOptExperiment = Model.optimizeCellCounts(fimResults,sum(Nc),'Determinant',[],[],10000*ones(size(Nc)));
            FIMOptExpt = Model.totalFim(fimResults,NcOptExperiment);
            Model.plotMHResults(MHResults,[FIM,FIMOptExpt])         

            % Find optimal experiment design given parameters sets but
            % where there is a base of 10 cells at every time point.
            NcBase = Nc;
            NcOptExperimentBase = Model.optimizeCellCounts(fimResults,sum(Nc),'Determinant',[],NcBase,10000*ones(size(Nc)));
            FIMOptExptBase = Model.totalFim(fimResults,NcOptExperimentBase+NcBase);
            Model.plotMHResults(MHResults,[FIM,FIMOptExpt,FIMOptExptBase])

        end 
        
        function LoadModelTest(testCase)
            % Test to see if the code can correctly load a model from file,
            % asscoiate data to that model, and then run a simple fitting
            % routine.
            %
            % Save Poisson Model to File
            delete('exampleResultsTest.mat')
            model = testCase.Poiss;
            save('TemporarySaveFile',"model")
            DataSettings = {'testData.csv',{'rna','exp1_s1'}};
            Pipeline = 'fittingPipelineExample';
            pipelineArgs.maxIter = 10;
            pipelineArgs.display = 'none';
            pipelineArgs.makePlot = true;
            saveFile = 'exampleResultsTest.mat';

            % Create model from preset, associate with data, run
            % 'fittingPipeline', and save result.
            SSIT('TemporarySaveFile','model',DataSettings,Pipeline,pipelineArgs,saveFile);

            % Load model from file, run 'fittingPipeline', and save result.
            SSIT(saveFile,'model',[],Pipeline,pipelineArgs,saveFile);

            testCase.verifyEqual(exist('exampleResultsTest.mat','file'), 2, ...
                'Model Creation Failed');

        end

        function TestAdvancedDataLoading(testCase)
            % Test loading multiple data sets and running logical data selection.
            model = testCase.Poiss;
            model = model.loadData({'test_data/fakeData4Testing1.xlsx','test_data/fakeData4Testing2.xlsx'},...
                {'rna','cyt'},...
                {'Replica',1,'>'});
            combineCellNum = model.dataSet.nCells;
            diff = max(abs(combineCellNum - [8;4]));
            testCase.verifyEqual(diff==0, true, ...
                'Advanced data loading resulted in incorrect cell number');

            model = model.loadData({'test_data/fakeData4Testing1.xlsx','test_data/fakeData4Testing2.xlsx'},...
                {'rna','cyt'},...
                {'Replica',1,'~='});
            combineCellNum = model.dataSet.nCells;
            diff = max(abs(combineCellNum - [8;4]));
            testCase.verifyEqual(diff==0, true, ...
                'Advanced data loading resulted in incorrect cell number');

            % Test custom constraint
            model = model.loadData({'test_data/fakeData4Testing1.xlsx','test_data/fakeData4Testing2.xlsx'},...
                {'rna','cyt'},...
                {[],[],'contains(TAB.Condition,''a'')&contains(TAB.Condition,''b'')'});
            combineCellNum = model.dataSet.nCells;
            diff = max(abs(combineCellNum - [2;1]));
            testCase.verifyEqual(diff==0, true, ...
                'Advanced data loading resulted in incorrect cell number');

            % Test manipulation of data (in this case summing two columns).
            model = model.loadData({'test_data/fakeData4Testing1.xlsx','test_data/fakeData4Testing2.xlsx'},...
                {'rna',[],'TAB.nuc+TAB.cyt'},...
                {[],[],'contains(TAB.Condition,''a'')&contains(TAB.Condition,''b'')'});
            distCounts = [model.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(1,4);
                model.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(1,5);
                model.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(1,6);];
            diff = max(abs(distCounts - [0;1;1]));
            testCase.verifyEqual(diff==0, true, ...
                'Datafield led to incorrec totals.');
            
        end
    end
end