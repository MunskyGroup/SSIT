classdef hybridModelTest < matlab.unittest.TestCase
    properties
         Hybrid
         HybridFSPSoln
    end
    methods (TestClassSetup)
        % Shared setup for the entire test class
         function setupModel(testCase)
            testCase.Hybrid = SSIT;  
            testCase.Hybrid.species = {'A';'B';'C'}; 
            testCase.Hybrid.initialCondition = [10;10;10];           
            testCase.Hybrid.propensityFunctions = {'k1*A';'k2*A';'k3*A*B';'k4*C'};         
            testCase.Hybrid.stoichiometry = [1,-1,0,0;
                                             0,1,-1,0;
                                             0,0,1,-1]; 
            testCase.Hybrid.parameters = ({'k1',0.6;'k2',0.2;'k3',0.3;'k4',0.4});  
            testCase.Hybrid.tSpan = linspace(0,2,21);
            testCase.Hybrid = testCase.Hybrid.formPropensitiesGeneral('Hybrid',true);
            testCase.Hybrid.summarizeModel;  
         end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function hybridModelWarningTest(testCase)
            testCase.Hybrid.useHybrid = true;
            testCase.Hybrid.hybridOptions.upstreamODEs = {'A','B'};
            testCase.Hybrid.summarizeModel
            % This should issue a warning and SSIT should ``automatically delete the upstream effect [...] from the stoichiometry for the downstream reaction"
            [testCase.HybridFSPSoln, testCase.Hybrid.fspOptions.bounds] = testCase.Hybrid.solve;
            testCase.Hybrid.summarizeModel
        end
    end
end