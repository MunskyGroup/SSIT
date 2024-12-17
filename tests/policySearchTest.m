classdef policySearchTest < burstyTest

    properties
         SensSoln
         fspFIM
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
         function setupModel(testCase)
             [testCase.Bursty, testCase.BurstyFSPSoln] = modelBuilder.buildBurstyModel();
             testCase.Bursty = testCase.Bursty.formPropensitiesGeneral('Bursty',true);
         end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)

        function PolicySearch(testCase)
            % Check that fminsearch resturns policy in log space.
            testCase.Bursty.solutionScheme = 'fspSens';
             testCase.Bursty.fspOptions.fspTol = 1e-6;
             OrigPars = [testCase.Bursty.parameters{:,2}]';
             for i = 1:3
                 testCase.Bursty.parameters(:,2) = num2cell(OrigPars.*(1+0.1*randn(size(OrigPars))));
                 testCase.SensSoln = testCase.Bursty.solve;
                 testCase.fspFIM = testCase.Bursty.computeFIM(testCase.SensSoln.sens);
             end

            BurstyPS = testCase.Bursty;
            BurstyPS.fspOptions.usePiecewiseFSP = true;
            BurstyPS.solutionScheme = 'FSP';
            BurstyPS.tSpan = [0,5,10,15];
            BurstyPS.propensityFunctions = {'kon*IGR*offGene';...
                'koff1*onGene*(t<5)+koff2*onGene*(t>=5)*(t<10)+koff3*onGene*(t>=15)';...
                'kr*onGene'; 'gr*rna'};
            Target = [0,50,50,50];
            BurstyPS.parameters = {'koff1', 0.02; 'koff2', 0.02; 'koff3', 0.02;...
                'kon', 0.002; 'kr', 1; 'gr', 0.005; 'a1', 2; 'r1', 0.05; 'r2', 0.1};
            logPolicy = zeros(1,3);

            loss = @(logPolicy)computeLossBursty(logPolicy,BurstyPS,Target);

            logPolicy = fminsearch(loss,logPolicy);

            fspPolicy = exp(logPolicy)
            [finalLoss,meanVal] = computeLossBursty(logPolicy,BurstyPS,Target);
            finalLoss

            plot(BurstyPS.tSpan,meanVal,BurstyPS.tSpan,Target,'x')

            function [loss,meanVal] = computeLossBursty(logPolicy,BurstyPS,Target)
                BurstyPS.parameters(1:3,2) = num2cell(exp(logPolicy));
                FSPsoln = BurstyPS.solve;
                meanVal = zeros(1,4);
                for i = 1:4
                    P = squeeze(sum(double(FSPsoln.fsp{i}.p.data),[1,2]));
                    meanVal(i) = (0:length(P)-1)*P;
                end
                loss = sum((Target - meanVal).^2);
            end
        end
    end
end