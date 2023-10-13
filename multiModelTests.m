classdef multiModelTests < matlab.unittest.TestCase
    % this tests:
    %       1) Creation of models using SBML.
    properties
        Poiss
        Poiss2
        MultiModel
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
            
            delete 'testData1.csv'
            tc.Poiss.ssaOptions.nSimsPerExpt = 100;
            tc.Poiss.ssaOptions.Nexp = 1;
            tc.Poiss.sampleDataFromFSP(tc.PoissSolution,'testData1.csv')
            tc.Poiss = tc.Poiss.loadData('testData1.csv',{'rna','exp1_s1'});

            delete 'testData2.csv'
            tc.Poiss2 = tc.Poiss;
            tc.Poiss2.ssaOptions.nSimsPerExpt = 100;
            tc.Poiss2.ssaOptions.Nexp = 1;
            tc.Poiss2.parameters = ({'kr',15;'gr',1});
            tc.Poiss2.sampleDataFromFSP(tc.PoissSolution,'testData2.csv')
            tc.Poiss2 = tc.Poiss2.loadData('testData2.csv',{'rna','exp1_s1'});

         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)

        function createMultiModel(tc)
            % Tests the loading of a model from SBML.
        end

    end
end