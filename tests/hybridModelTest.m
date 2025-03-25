classdef hybridModelTest < matlab.unittest.TestCase
    properties
         Bursty
         BurstyFSPSoln
    end
    methods (TestClassSetup)
        % Shared setup for the entire test class
         function setupModel(testCase)
            testCase.Bursty = SSIT;  
            testCase.Bursty.species = {'offGene';'onGene';'rna'}; 
            testCase.Bursty.initialCondition = [1;0;0];           
            testCase.Bursty.propensityFunctions = {'kon*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
            testCase.Bursty.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'}; 
            testCase.Bursty.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; 
            testCase.Bursty.parameters = ({'kon',0.6;'koff',2;'kr',0.3;'gr',0.04});  
            testCase.Bursty.summarizeModel;  
         end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function hybridModelWarningTest(testCase)
            testCase.Bursty.useHybrid = true;
            testCase.Bursty.hybridOptions.upstreamODEs = {'offGene','onGene'};
            testCase.Bursty.summarizeModel
            [testCase.BurstyFSPSoln, testCase.Bursty.fspOptions.bounds] = testCase.Bursty.solve;
        end
    end
end