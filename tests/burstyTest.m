classdef burstyTest < matlab.unittest.TestCase
    properties
        Bursty
        BurstyFSPSoln
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(testCase1)
            addpath(genpath('../src'));
            testCase1.Bursty = SSIT;  
            testCase1.Bursty.species = {'offGene';'onGene';'rna'}; 
            testCase1.Bursty.initialCondition = [2;0;0];           
            testCase1.Bursty.propensityFunctions = {'kon*IGR*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
            testCase1.Bursty.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'}; 
            testCase1.Bursty.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; 
            testCase1.Bursty.parameters = ({'koff',0.014;'kon',0.002;'kr',1;'gr',0.004;...
                                            'a1',20;'r1',0.04;'r2',0.1});  
            testCase1.Bursty.fspOptions.initApproxSS = true;  % Set Initial Distribution to Steady State.
            testCase1.Bursty.summarizeModel                   % Print visual summary of model            
            [testCase1.BurstyFSPSoln,testCase1.Bursty.fspOptions.bounds] = testCase1.Bursty.solve;
            tic
            [testCase1.BurstyFSPSoln,testCase1.Bursty.fspOptions.bounds] = testCase1.Bursty.solve(testCase1.BurstyFSPSoln.stateSpace); 
            testCase1.BurstyFSPSoln.time = toc;

            delete 'testData.csv'
            testCase1.Bursty.ssaOptions.nSimsPerExpt = 1000;
            testCase1.Bursty.ssaOptions.Nexp = 1;
            testCase1.Bursty.sampleDataFromFSP(testCase1.BurstyFSPSoln,'testData.csv');

            testCase1.Bursty = testCase1.Bursty.loadData('testData.csv',{'rna','exp1_s1'});
         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function ModelCreation(testCase)
            % In this trivial test, we check that the SSIT is set up with
            % the right names for the 'rna' species.
            nm = testCase.Bursty.species;
            testCase.verifyEqual(nm{1}, 'offGene', ...
                'Species name is incorrect');
            testCase.verifyEqual(nm{2}, 'onGene', ...
                'Species name is incorrect');
            testCase.verifyEqual(nm{3}, 'rna', ...
                'Species name is incorrect');
        end

        function FspConverged(testCase)
            % In this test, we check tha tthe FSP solution exits with an
            % appropriate FSP tolerance value.
            final = testCase.BurstyFSPSoln.fsp{end}.p.sum;
            tst = (1-final)<=testCase.Bursty.fspOptions.fspTol;
            testCase.verifyEqual(tst, true, ...
                'Final FSP is not within tolerance');
        end
    end
end