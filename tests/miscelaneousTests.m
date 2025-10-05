classdef miscelaneousTests < matlab.unittest.TestCase
    % this tests:
    %       1) Creation of models using SBML.
    properties
        Model
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath(genpath('../src'));
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
            tc.Model = tc.Model.formPropensitiesGeneral('SBMEModel');
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
            TwoDNonLinearTV = TwoDNonLinearTV.formPropensitiesGeneral('TwoNonLinTV',false);
            [fspSoln1] = TwoDNonLinearTV.solve;
            TwoDNonLinearTV.makePlot(fspSoln1,'meansAndDevs',[],[],[1])
            
            figure(2)
            TwoDNonLinearTV.exportSimBiol(true);
            TwoDNonLinearTV.exportToSBML('TwoDNonLinearTV.xml');

            clear TwoDNonLinearTV

            TwoDNonLinearTV = SSIT;
            TwoDNonLinearTV = TwoDNonLinearTV.createModelFromSBML('TwoDNonLinearTV.xml');
            TwoDNonLinearTV.tSpan = linspace(0,10,40);
                        
            TwoDNonLinearTV = TwoDNonLinearTV.formPropensitiesGeneral('TwoNonLinTV',false);
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
            FourDNonLinearTV = SSIT;
            FourDNonLinearTV.species = {'rna','prot','phosProt','GFP'};
            FourDNonLinearTV.initialCondition = [0;0;0;0];
            FourDNonLinearTV.propensityFunctions = {'kr';'gr*rna*(1/(1+(phosProt/M)^eta*Ir))';...
                'kp*rna';'gp*prot';...
                'kx*prot';'gx*phosProt';...
                'kg*rna';'ggfp*GFP'};
            FourDNonLinearTV.stoichiometry = [1,-1,0,0,0,0,0,0;...
                0,0,1,-1,-1,0,0,0;...
                0,0,0,0,1,-1,0,0;...
                0,0,0,0,0,0,1,-1];
            FourDNonLinearTV.parameters = ({'kr',20;'gr',1;'kp',18;'gp',1;'M',20;'eta',5;...
                'kx',2;'gx',1;'kg',12;'ggfp',0.5});
            FourDNonLinearTV.inputExpressions = {'Ir','1+cos(t)'};
            FourDNonLinearTV.tSpan = linspace(0,10,40);
            FourDNonLinearTV.solutionScheme = 'ode';
            FourDNonLinearTV = FourDNonLinearTV.formPropensitiesGeneral('TwoDTV',false);
            
            FourDNonLinearTV.odeIntegrator = 'ode23s';

            [odeSoln1] = FourDNonLinearTV.solve;
            parVector = [FourDNonLinearTV.parameters{:,2}];

            Ntests = 100;
            parVectorSets = repmat(parVector,Ntests,1).*(1+0.1*randn(Ntests,size(parVector,2)));
            results = zeros(Ntests,4);
            resultsSB = zeros(Ntests,4);

            tic
            for i=1:Ntests
                FourDNonLinearTV.parameters(:,2) = num2cell(parVectorSets(i,:));
                [odeSoln1] = FourDNonLinearTV.solve;
                results(i,:) = odeSoln1.ode(end,:);
            end
            SSITSolveTime100pars = toc

            sbModel = FourDNonLinearTV.exportSimBiol;
            csObj = getconfigset(sbModel,'active');
            set(csObj,'Stoptime',max(FourDNonLinearTV.tSpan));
            
            [t,x,names] = sbiosimulate(sbModel);
            tic
            for i=1:Ntests
                for j=1:length(parVector)
                    sbModel.Parameters(j).Value = parVectorSets(i,j);
                end
                [t,x,names] = sbiosimulate(sbModel);
                resultsSB(i,:) = x(end,:);
            end
            simBiolSolveTime100pars = toc

            tc.verifyEqual(SSITSolveTime100pars<(2*simBiolSolveTime100pars), true, ...
                'SSIT ODE Solution is > 2x slower than SimBiology.');

            meanError = mean(abs((resultsSB-results)./(mean(results))),"all");
            
            tc.verifyEqual(meanError<0.01, true, ...
                'Average SSIT ODE Solution is not within 1 percent of SimBiology.');

            % TODO -- Need to add SSIT option for ODE aanlysis for user to
            % choose the ODE solver.  For stiff problems, ODE23s is much
            % faster, but for non-stiff problems (like this one) ODE45 is
            % faster. 

         end

         function testMoments(tc)
             close all
             ThreeDNLinearTV = SSIT;
             ThreeDNLinearTV.species = {'rna','prot','phosProt'};
             ThreeDNLinearTV.initialCondition = [0;0;0];
             ThreeDNLinearTV.propensityFunctions = {'kr';'gr*rna*Ir';...
                 'kp*rna';'gp*prot';...
                 'kx*prot';'gx*phosProt'};
             ThreeDNLinearTV.stoichiometry = [1,-1,0,0,0,0;...
                 0,0,1,-1,-1,0;...
                 0,0,0,0,1,-1];
             ThreeDNLinearTV.parameters = ({'kr',5;'gr',1;'kp',3;'gp',1;'M',20;'eta',5;...
                 'kx',2;'gx',1});
             ThreeDNLinearTV.inputExpressions = {'Ir','1+cos(t)'};
             ThreeDNLinearTV.tSpan = linspace(0,10,40);

             ThreeDNLinearTV.odeIntegrator = 'ode45';

             ThreeDNLinearTV.solutionScheme = 'ode';
             ThreeDNLinearTV = ThreeDNLinearTV.formPropensitiesGeneral('TwoDTV',false);
             [~,~,ThreeDNLinearTV] = ThreeDNLinearTV.solve;

             ThreeDNLinearTV.solutionScheme = 'moments';
             ThreeDNLinearTV = ThreeDNLinearTV.formPropensitiesGeneral('TwoDTV',false);
             [~,~,ThreeDNLinearTV] = ThreeDNLinearTV.solve;

             maxRelError = max((ThreeDNLinearTV.Solutions.moments(1:3,:) - ThreeDNLinearTV.Solutions.ode')./ThreeDNLinearTV.Solutions.ode',[],'all');

             tc.verifyEqual(maxRelError<0.01, true, ...
                'ODE and moments mean solutions do not match.');


             ThreeDNLinearTV.solutionScheme = 'fsp';
             ThreeDNLinearTV.fspOptions.verbose = true;
             ThreeDNLinearTV = ThreeDNLinearTV.formPropensitiesGeneral('TwoDTV',false);
             [~,~,ThreeDNLinearTV] = ThreeDNLinearTV.solve;

             % sum(double(ThreeDNLinearTV.Solutions.fsp{end}.p.data),[2,3])ThreeDNLinearTV.Solutions.stateSpace.states
             % finalMean = ThreeDNLinearTV.Solutions.fsp{end}.p

             for i = 1:length(ThreeDNLinearTV.Solutions.fsp)
                 sz = size(ThreeDNLinearTV.Solutions.fsp{i}.p.data);
                 fspMn1(i) = [0:sz(1)-1]*double(ThreeDNLinearTV.Solutions.fsp{i}.p.sumOver([2,3]).data);
                 fspMn2(i) = [0:sz(2)-1]*double(ThreeDNLinearTV.Solutions.fsp{i}.p.sumOver([1,3]).data);
                 fspMn3(i) = [0:sz(3)-1]*double(ThreeDNLinearTV.Solutions.fsp{i}.p.sumOver([1,2]).data);
                 fspMn1_2(i) = [0:sz(1)-1].^2*double(ThreeDNLinearTV.Solutions.fsp{i}.p.sumOver([2,3]).data);
                 fspMn2_2(i) = [0:sz(2)-1].^2*double(ThreeDNLinearTV.Solutions.fsp{i}.p.sumOver([1,3]).data);
                 fspMn3_2(i) = [0:sz(3)-1].^2*double(ThreeDNLinearTV.Solutions.fsp{i}.p.sumOver([1,2]).data);
             end

             % plot(fspMn3_2); hold on;
             % plot(ThreeDNLinearTV.Solutions.moments(end,:),'x')

             relErr = max((fspMn3_2 - ThreeDNLinearTV.Solutions.moments(end,:))./fspMn3_2,[],'all');

             tc.verifyEqual(relErr<0.01, true, ...
                 'Relative error of second moment is not below 1%.');



         end
         
         function testCommandLine(testCase)
             %% Define Model/Data Combination
             % Specify data set to fit.
             DataSettings = {'testData.csv',{'mRNA','exp1_s1'}};

             % Create model from preset, associate with data.
             Model = SSIT('BirthDeath',[],DataSettings);

             % Save model for later use.
             saveFile = 'temporarySaveFile';
             save(saveFile,"Model")

             %% Call Pipeline to Fit Model
             % Specify pipeline to apply to model and arguments
             % ("../src/exampleData/examplePipelines/fittingPipelineExample.m")
             Pipeline = 'fittingPipelineExample';
             pipelineArgs.maxIter = 10;
             pipelineArgs.display = 'iter';
             pipelineArgs.makePlot = false;

             saveFile2 = 'temporarySaveFile2.mat';

             if exist(saveFile2,"file")
                 delete(saveFile2);
             end
             if exist('logFile.txt','file')
                 delete('logFile.txt');
             end

             SSIT.generateCommandLinePipeline(saveFile,'Model',[],Pipeline,pipelineArgs,saveFile2,'logFile.txt',true)

             disp('Pausing to give time for background run to complete.')
             tic
             success = false;
             while toc<90
                 if exist(saveFile2,'file')&&exist('logFile.txt','file')
                     disp(['Background run completed in ',num2str(toc),' s'])
                     success =true;
                     break
                 end
             end
             testCase.verifyEqual(success, true, ...
                 'Files from background fit were not created.');

         end

         

    end
end