classdef modelBuilder
    methods (Static)
        function [Bursty, BurstyFSPSoln] = buildBurstyModel()
            addpath(genpath('../src'));
            Bursty = SSIT;  
            Bursty.species = {'offGene';'onGene';'rna'}; 
            Bursty.initialCondition = [1;0;0];           
            Bursty.propensityFunctions = {'kon*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
            %Bursty.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'}; 
            Bursty.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; 
            Bursty.parameters = ({'kon',0.6;'koff',2;'kr',0.3;'gr',0.04});  
            Bursty.summarizeModel;                         
            [BurstyFSPSoln,Bursty.fspOptions.bounds] = Bursty.solve;
            [BurstyFSPSoln,Bursty.fspOptions.bounds] = Bursty.solve(BurstyFSPSoln.stateSpace); 
        end
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

            [testCase1.PoissSolution,testCase1.Poiss.fspOptions.bounds] = testCase1.Poiss.solve;
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
end