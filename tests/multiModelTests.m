classdef multiModelTests < matlab.unittest.TestCase
    % this tests:
    %       1) Creation of models using SBML.
    properties
        Poiss
        PoissSolution
        Poiss2
        Poiss2Solution
        MultiModel
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestModel1(tc)
            addpath(genpath('../src'));
            %% Test Case 1 - a simple Poisson model
            tc.Poiss = SSIT;
            tc.Poiss.species = {'rna'};
            tc.Poiss.initialCondition = 0;
            tc.Poiss.propensityFunctions = {'kr';'gr*rna'};
            tc.Poiss.stoichiometry = [1,-1];
            tc.Poiss.parameters = ({'kr',10;'gr',1});
            tc.Poiss.tSpan = linspace(0,2,9);
            tc.Poiss.fspOptions.fspTol = 1e-5;
            tc.Poiss = tc.Poiss.formPropensitiesGeneral('Poiss1',true);
            tc.Poiss.fittingOptions.modelVarsToFit = [1,2];
            [tc.PoissSolution,tc.Poiss.fspOptions.bounds] = tc.Poiss.solve;
            
            tc.Poiss.ssaOptions.nSimsPerExpt = 100;
            tc.Poiss.ssaOptions.Nexp = 1;
            delete 'testData1.csv'
            tc.Poiss.sampleDataFromFSP(tc.PoissSolution,'testData1.csv');
            tc.Poiss = tc.Poiss.loadData('testData1.csv',{'rna','exp1_s1'});

            tc.Poiss2 = tc.Poiss;
            tc.Poiss2.ssaOptions.nSimsPerExpt = 100;
            tc.Poiss2.ssaOptions.Nexp = 1;
            tc.Poiss2.parameters = ({'kr',15;'gr',1});
            tc.Poiss2 = tc.Poiss2.formPropensitiesGeneral('Poiss2',true);
            [tc.Poiss2Solution,tc.Poiss2.fspOptions.bounds] = tc.Poiss2.solve;
            delete 'testData2.csv'
            tc.Poiss2.sampleDataFromFSP(tc.Poiss2Solution,'testData2.csv');
            tc.Poiss2 = tc.Poiss2.loadData('testData2.csv',{'rna','exp1_s1'});

         end  
    end

    methods (TestMethodSetup)
        % Setup for each test

    end

    methods (Test)

        function createMultiModel(tc)
            % Tests the creation, fitting, and saving of a combined model
            parset = {[1,2],[3,2]};
            ModelSet = {tc.Poiss,tc.Poiss2};
            combinedModel = SSITMultiModel(ModelSet,parset);
            combinedModel = combinedModel.initializeStateSpaces;
            parGuess = [1,1,1];

            % Fit the model to the data
            fitOptions = optimset('Display','final','MaxIter',500);
            parGuess = combinedModel.maximizeLikelihood(...
                parGuess, fitOptions);

            combinedModel = combinedModel.updateModels(parGuess,true);

            delete('temporaryHybridModel.mat')
            save('temporaryHybridModel',"combinedModel")
        end

        function loadHybridModelFromFile(tc)
            delete('exampleResultsTest.mat')
            
            DataSettings = {'testData1.csv',{'rna','exp1_s1'},{};...
                'testData2.csv',{'rna','exp1_s1'},{}};
            
            Pipeline = 'multiModelFittingPipelineExample';
            pipelineArgs.maxIter = 10;
            pipelineArgs.display = 'iter';
            saveFile = 'exampleCombinedResultsTest.mat';
            delete(saveFile)
            
            SSIT('temporaryHybridModel','combinedModel',DataSettings,Pipeline,pipelineArgs,saveFile);
            
            tc.verifyEqual(exist(saveFile,'file'), 2, ...
                'Combined Model Creation Failed');        
        end

        function testFIM(tc)
            parset = {[1,2],[3,2]};
            ModelSet = {tc.Poiss,tc.Poiss2};
            combinedModel = SSITMultiModel(ModelSet,parset);
            combinedModel = combinedModel.initializeStateSpaces;
            parGuess = [tc.Poiss.parameters{:,2},tc.Poiss2.parameters{1,2}];

            combinedModel = computeFIMs(combinedModel);
            
            t = tc.Poiss.tSpan;
            k1 = tc.Poiss.parameters{1,2};
            k2 = tc.Poiss2.parameters{1,2};
            g = tc.Poiss.parameters{2,2};
            
            lam1 = k1/g*(1-exp(-g*t))';
            lam2 = k2/g*(1-exp(-g*t))';

            exactTotalFIM = zeros(3);
            
            for i = 2:length(t)
                exactIk1(i) = lam1(i)/k1^2;
                exactIk2(i) = lam2(i)/k2^2;
                
                fspIk1(i) = combinedModel.FIM.fims{1}{i}(1,1);
                fspIk2(i) = combinedModel.FIM.fims{2}{i}(1,1);
                
                exactIg1(i) = (-lam1(i)/g + k1*t(i)/g*exp(-g*t(i)))^2/lam1(i);
                exactIg2(i) = (-lam2(i)/g + k2*t(i)/g*exp(-g*t(i)))^2/lam2(i);
                
                fspIg1(i) = combinedModel.FIM.fims{1}{i}(2,2);
                fspIg2(i) = combinedModel.FIM.fims{2}{i}(2,2);
                
                exactIkg1(i) = (-lam1(i)/g + k1*t(i)/g*exp(-g*t(i)))/k1;
                exactIkg2(i) = (-lam2(i)/g + k2*t(i)/g*exp(-g*t(i)))/k2;
                
                fspIkg1(i) = combinedModel.FIM.fims{1}{i}(1,2);
                fspIkg2(i) = combinedModel.FIM.fims{2}{i}(1,2);

                exactTotalFIM = exactTotalFIM + tc.Poiss.ssaOptions.nSimsPerExpt*...
                    [exactIk1(i),exactIkg1(i),0;...
                    exactIkg1(i),exactIg1(i)+exactIg2(i),exactIkg2(i);...
                    0,exactIkg2(i),exactIk2(i)];
            end
            % 
            diff = max(abs((exactIk1-fspIk1)./exactIk1)+...
                abs((exactIg1-fspIg1)./exactIg1)+...
                abs((exactIkg1-fspIkg1)./exactIkg1))+...
                max(abs((exactIk2-fspIk2)./exactIk2)+...
                abs((exactIg2-fspIg2)./exactIg2)+...
                abs((exactIkg2-fspIkg2)./exactIkg2));
           
            diff = max(diff,max(abs((combinedModel.FIM.totalFIM - exactTotalFIM)./exactTotalFIM),[],"all")/9);
            
            tc.verifyEqual(diff<0.001, true, ...
                'FIM Calculation is not within 0.1% Tolerance');
        end

    end
end