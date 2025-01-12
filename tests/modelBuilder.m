classdef modelBuilder
    methods (Static)
        function [Bursty, BurstyFSPSoln] = buildBurstyModel()
            addpath(genpath('../src'));
            %% a simple Bursting Gene model
            Bursty = SSIT;  
            Bursty.species = {'offGene';'onGene';'rna'}; 
            Bursty.initialCondition = [1;0;0];           
            Bursty.propensityFunctions = {'kon*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
            Bursty.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'}; 
            Bursty.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; 
            Bursty.parameters = ({'kon',0.6;'koff',2;'kr',0.3;'gr',0.04});  
            Bursty.summarizeModel;                         
            [BurstyFSPSoln,Bursty.fspOptions.bounds] = Bursty.solve;
            [BurstyFSPSoln,Bursty.fspOptions.bounds] = Bursty.solve(BurstyFSPSoln.stateSpace); 
        end
        function [Poiss, PoissSolution, PoissODE] = buildPoissModel()
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
        function [Tog, TogFSPSoln, TogfspFIM, TogSensSoln] = buildTogSwitchModel()
            addpath(genpath('../src'));
            %% a simple Toggle Switch model
            Tog = SSIT;
            Tog.species = {'lacI'; 'lambdaCI'};
            Tog.initialCondition = [4;0];
            Tog.parameters = {'p1', 0.7455; 'p2', 0.3351; 'p3', 0.0078; 'p4', 0.4656; 'p5', 0.0193; ...
                                'p6', 0.2696; 'p7', 2.5266; 'p8', 0.4108; 'p9', 0.6880; ...
                                'xi1', 0.9276; 'xi2', 0.2132; 'xi3', 0.8062; 'xi4', 0.3897};
            Tog.inputExpressions = {'IUV','xi1*(t<5)+xi2*(t>=5)*(t<10)+xi3*(t>=10)*(t<15)+xi4*(t>=15)' };
            Tog.fspOptions.usePiecewiseFSP = true;
            Tog.propensityFunctions = {'p1+(p3/(1+(p6*(lambdaCI)^p8)))';
                                    'p2+(p4/(1+(p5*(lacI)^p7)))'; 'p9*lacI'; 'IUV*lambdaCI'};
            Tog.stoichiometry = [1,0,-1,0; 0,1,0,-1];
            Tog.tSpan = [0,5,10,15,20];
            Tog.fspOptions.fspTol = 1e-5;
            Tog = Tog.formPropensitiesGeneral('Tog', true);
            Tog.solutionScheme = 'FSP';
            Tog.fspOptions.fspTol = 1e-7;
            [TogFSPSoln, Tog.fspOptions.bounds] = Tog.solve;
            Tog.solutionScheme = 'fspSens';
            Tog.fspOptions.fspTol = 1e-6;
            OrigPars = [Tog.parameters{:,2}]';
            Tog.sensOptions.solutionMethod = 'finiteDifference';
            for i = 1:10
                Tog.parameters(:,2) = num2cell(OrigPars.*(1+0.1*randn(size(OrigPars))));
                [TogSensSoln, Tog.fspOptions.bounds] = Tog.solve;
                TogfspFIM = Tog.computeFIM(SensSoln.sens);
            end
        end
    end
end