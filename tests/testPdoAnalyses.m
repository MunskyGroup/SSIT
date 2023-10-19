classdef testPdoAnalyses < matlab.unittest.TestCase
    properties
        Poiss
        ModelPDOSpots
        Solution
        SensSoln
        MetHastResults
        
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath('../CommandLine')
            %% Test Case 1 - a simple Poisson model
            tc.Poiss = SSIT;
            tc.Poiss.species = {'rna'};
            tc.Poiss.initialCondition = 0;
            tc.Poiss.propensityFunctions = {'kr';'gr*rna'};
            tc.Poiss.stoichiometry = [1,-1];
            tc.Poiss.parameters = ({'kr',10;'gr',1});
            tc.Poiss.tSpan = linspace(0,2,21);
            tc.Poiss.fspOptions.fspTol = 1e-5;
            tc.Poiss = tc.Poiss.formPropensitiesGeneral('Poiss');

            %% Compute CME Solution
            tc.Poiss.solutionScheme = 'FSP';
            tc.Poiss.fspOptions.fspTol = 1e-6;            
            [tc.Solution.UnCorrected,tc.Poiss.fspOptions.bounds] = tc.Poiss.solve;

            %% Calibrate on Empirical Data to Find PDO
            tc.ModelPDOSpots = tc.Poiss.calibratePDO('../ExampleData/pdoCalibrationData.csv',...
                {'rna'},{'nTotal'},{'nSpots0'},'AffinePoiss',true);

            %% Generate Distorted Data
            tc.Solution.Corrected = tc.ModelPDOSpots.solve;
            delete 'testDataDistorted.csv'
            tc.ModelPDOSpots.ssaOptions.nSimsPerExpt = 1000;
            tc.ModelPDOSpots.ssaOptions.Nexp = 1;
            tc.ModelPDOSpots.sampleDataFromFSP(tc.Solution.Corrected,'testDataDistorted.csv')

            %% Compute CME Sensitivity
            tc.Poiss.solutionScheme = 'fspSens';
            tc.Poiss.fspOptions.fspTol = inf;            
            tc.SensSoln = tc.Poiss.solve;

            %% Associate Model with Data
            tc.Poiss = tc.Poiss.loadData('testDataDistorted.csv',{'rna','exp1_s1_Distorted'});
            tc.ModelPDOSpots = tc.ModelPDOSpots.loadData('testDataDistorted.csv',{'rna','exp1_s1_Distorted'});

            % tc.Poiss.computeLikelihood
            % tc.ModelPDOSpots.computeLikelihood
            
            %% Run MetHastings to find parameter uncertainty with and without distortion
            tc.Poiss.solutionScheme = 'FSP';
            tc.ModelPDOSpots.solutionScheme = 'FSP';

            % run Metropolis Hastings
            MHFitOptions.thin=1;
            MHFitOptions.numberOfSamples=1000;
            MHFitOptions.burnIn=0;
            MHFitOptions.progress=true;
            MHFitOptions.useFIMforMetHast = true;
            MHFitOptions.CovFIMscale = 1.0;
            MHFitOptions.numChains = 1;
            MHFitOptions.saveFile = 'TMPMHChain.mat';
            
            delete('TMPMHChain.mat')
            tc.Poiss.fittingOptions.modelVarsToFit = [1,2];
            [~,~,tc.MetHastResults.UnCorrected] = tc.Poiss.maximizeLikelihood(...
                [], MHFitOptions, 'MetropolisHastings');
            delete('TMPMHChain.mat')

            tc.ModelPDOSpots.fittingOptions.modelVarsToFit = [1,2];
            [~,~,tc.MetHastResults.Corrected] = tc.ModelPDOSpots.maximizeLikelihood(...
                [], MHFitOptions, 'MetropolisHastings');
            delete('TMPMHChain.mat')

        end
    end
    methods (Test)
        % Test methods

        function tc = ModelCreation(tc)
            % In this trivial test, we check that the SSIT is set up with
            % the right names for the 'rna' species.
            nm = tc.Poiss.species;
            tc.verifyEqual(nm{1}, 'rna', ...
                'Species name is incorrect');
           
        end

        function tc = FitModelWithPDO(tc)
        end

        function MetHastWithPDO(tc)
        end

        function SensitivityWithPDO(tc)
        end

        function FIMwithPDO(tc)
            fimsPDOSpot = tc.ModelPDOSpots.computeFIM(sensSoln);
            fimPDOSpots = tc.ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt);
            Model.plotMHResults(mhResults,{FIM,fimOpt,fimPDOSpots});
        end

    end

end