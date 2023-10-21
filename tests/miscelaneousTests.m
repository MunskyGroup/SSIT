classdef miscelaneousTests < matlab.unittest.TestCase
    % this tests:
    %       1) Creation of models using SBML.
    properties
        Model
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath('../CommandLine')
            tc.Model = SSIT;
         end  
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)

        function loadModelFromSBML(tc)
            % Tests the loading of a model from SBML.
            tc.Model = tc.Model.createModelFromSBML('../SBML_test_cases/00010/00010-sbml-l1v2.xml',true);
            [fspSoln] = tc.Model.solve;
            tc.Model.makePlot(fspSoln,'meansAndDevs')
        end

        function testNonlinearLogicTimeVarying(tc)

                        %% Test Case 3 - 2 Species Poisson Model
            TwoDNonLinearTV = SSIT;
            TwoDNonLinearTV.species = {'rna1','rna2'};
            TwoDNonLinearTV.initialCondition = [0;0];
            TwoDNonLinearTV.propensityFunctions = {'kr1';'gr1*rna1*(1/(1+(rna2/M)^eta*Ir))';'kr2';'gr2*rna2'};
            TwoDNonLinearTV.stoichiometry = [1,-1,0,0;0,0,1,-1];
            TwoDNonLinearTV.parameters = ({'kr1',20;'gr1',1;'kr2',18;'gr2',1;'M',20;'eta',5});
            TwoDNonLinearTV.inputExpressions = {'Ir','t>1'};
            TwoDNonLinearTV.tSpan = linspace(0,2,21);
            TwoDNonLinearTV = TwoDNonLinearTV.formPropensitiesGeneral('TwoNonLinTV');

            [TwoDNonLinearTVSolution,TwoDNonLinearTV.fspOptions.bounds] = TwoDNonLinearTV.solve;
            tic
            [TwoDNonLinearTVSolution,TwoDNonLinearTV.fspOptions.bounds] = TwoDNonLinearTV.solve(TwoDNonLinearTVSolution.stateSpace);
            TwoDNonLinearTVSolution.time = toc;

            % %% ODE model for Two Poisson process
            % TwoDNonLinearTVODE = TwoDNonLinearTV;
            % TwoDNonLinearTVODE.solutionScheme = 'ode';
            % TwoDNonLinearTVODE = TwoDNonLinearTVODE.formPropensitiesGeneral('TwoDNonLinearTVODE');  


        end

        function ExportToSBML(tc)
            close all
            TwoDNonLinearTV = SSIT;
            TwoDNonLinearTV.species = {'rna1','rna2'};
            TwoDNonLinearTV.initialCondition = [0;0];
            TwoDNonLinearTV.propensityFunctions = {'kr1';'gr1*rna1*(1/(1+(rna2/M)^eta*Ir))';'kr2';'gr2*rna2'};
            TwoDNonLinearTV.stoichiometry = [1,-1,0,0;0,0,1,-1];
            TwoDNonLinearTV.parameters = ({'kr1',20;'gr1',1;'kr2',18;'gr2',1;'M',20;'eta',5});
            TwoDNonLinearTV.inputExpressions = {'Ir','1+cos(t)'};
            TwoDNonLinearTV.tSpan = linspace(0,10,40);
            [fspSoln1] = TwoDNonLinearTV.solve;
            TwoDNonLinearTV.makePlot(fspSoln1,'meansAndDevs',[],[],[1])
            
            figure(2)
            TwoDNonLinearTV.exportSimBiol(true);
            TwoDNonLinearTV.exportToSBML('TwoDNonLinearTV.xml');

            clear TwoDNonLinearTV

            TwoDNonLinearTV = SSIT;
            TwoDNonLinearTV = TwoDNonLinearTV.createModelFromSBML('TwoDNonLinearTV.xml');
            TwoDNonLinearTV.tSpan = linspace(0,10,40);
            [fspSoln2] = TwoDNonLinearTV.solve;
            TwoDNonLinearTV.makePlot(fspSoln2,'meansAndDevs',[],[],[3])

            P1 = double(fspSoln1.fsp{end}.p.data);
            P2 = double(fspSoln2.fsp{end}.p.data);
            diff = sum(abs(P1-P2),"all");
            tst = diff<1e-6;

            tc.verifyEqual(tst, true, ...
                'Exporting and Importing to SBML Changed final FSP solution');
            
        end

         function test2DSpeedVsSimBiol(tc)
            close all
            TwoDNonLinearTV = SSIT;
            TwoDNonLinearTV.species = {'rna','prot','phosProt','GFP'};
            TwoDNonLinearTV.initialCondition = [0;0;0;0];
            TwoDNonLinearTV.propensityFunctions = {'kr';'gr*rna*(1/(1+(phosProt/M)^eta*Ir))';...
                'kp*rna';'gp*prot';...
                'kx*prot';'gx*phosProt';...
                'kg*rna';'ggfp*GFP'};
            TwoDNonLinearTV.stoichiometry = [1,-1,0,0,0,0,0,0;...
                0,0,1,-1,-1,0,0,0;...
                0,0,0,0,1,-1,0,0;...
                0,0,0,0,0,0,1,-1];
            TwoDNonLinearTV.parameters = ({'kr',20;'gr',1;'kp',18;'gp',1;'M',20;'eta',5;...
                'kx',2;'gx',1;'kg',12;'ggfp',0.5});
            TwoDNonLinearTV.inputExpressions = {'Ir','1+cos(t)'};
            TwoDNonLinearTV.tSpan = linspace(0,10,40);
            TwoDNonLinearTV.solutionScheme = 'ode';
            TwoDNonLinearTV = TwoDNonLinearTV.formPropensitiesGeneral('TwoDTV');
            
            [odeSoln1] = TwoDNonLinearTV.solve;
            parVector = [TwoDNonLinearTV.parameters{:,2}];
            tic
            for i=1:100
                TwoDNonLinearTV.parameters(:,2) = num2cell(parVector.*(1+0.1*randn(size(parVector))));
                [odeSoln1] = TwoDNonLinearTV.solve;
            end
            SSITSolveTime100pars = toc

            sbModel = TwoDNonLinearTV.exportSimBiol;
            csObj = getconfigset(sbModel,'active');
            set(csObj,'Stoptime',max(TwoDNonLinearTV.tSpan));
            
            [t,x,names] = sbiosimulate(sbModel);
            tic
            for i=1:100
                for j=1:length(parVector)
                    sbModel.Parameters(j).Value = parVector(j)*(1+0.1*randn);
                end
                [t,x,names] = sbiosimulate(sbModel);
            end
            simBiolSolveTime100pars = toc

            tc.verifyEqual(SSITSolveTime100pars<(2*simBiolSolveTime100pars), true, ...
                'SSIT ODE Solution is > 2x slower than SimBiology.');
            
         end

    end
end