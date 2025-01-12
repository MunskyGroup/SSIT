classdef burstyTest < matlab.unittest.TestCase
    properties
         Bursty
         BurstyFSPSoln
    end
    methods (TestClassSetup)
        % Shared setup for the entire test class
         function setupModel(testCase)
             [testCase.Bursty, testCase.BurstyFSPSoln] = modelBuilder.buildBurstyModel();
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