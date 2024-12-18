classdef modelBuilder
    methods (Static)
        function [Bursty, BurstyFSPSoln] = buildBurstyModel()
            addpath(genpath('../src'));
            %% a simple Bursting Gene model
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
        function [Poiss, PoissSolution] = buildPoissModel()
            addpath(genpath('../src'));
            %% a simple Poisson (Birth-Death) model
            Poiss = SSIT;
            Poiss.species = {'rna'};
            Poiss.initialCondition = 0;
            Poiss.propensityFunctions = {'kr';'gr*rna'};
            Poiss.stoichiometry = [1,-1];
            Poiss.parameters = ({'kr',10;'gr',1});
            Poiss.tSpan = linspace(0,2,21);
            Poiss.fspOptions.fspTol = 1e-5;
            Poiss = Poiss.formPropensitiesGeneral('Poiss',true);

            [PoissSolution,Poiss.fspOptions.bounds] = Poiss.solve;
            tic
            [PoissSolution,Poiss.fspOptions.bounds] = Poiss.solve(PoissSolution.stateSpace);
            PoissSolution.time = toc;

            delete 'testData.csv'
            Poiss.ssaOptions.nSimsPerExpt = 1000;
            Poiss.ssaOptions.Nexp = 1;
            Poiss.sampleDataFromFSP(PoissSolution,'testData.csv');

            Poiss = Poiss.loadData('testData.csv',{'rna','exp1_s1'});

            %% ODE model for Poisson process
            PoissODE = Poiss;
            PoissODE.solutionScheme = 'ode';
            PoissODE = PoissODE.formPropensitiesGeneral('PoissODE');  
         end  
    end
end