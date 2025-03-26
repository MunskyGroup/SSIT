classdef pdoTest < burstyTest
    properties
        ModelPDO_Intensity
        ModelPDO_single_Z
        ModelPDO_5thru9_Zs
    end

    methods (TestClassSetup)
        % Define SSIT model	
        function setupModel(testCase)
            [testCase.Bursty, testCase.BurstyFSPSoln] = modelBuilder.buildBurstyModel();
            testCase.Bursty = testCase.Bursty.loadData('test_data/pdo.csv',{'rna','RNA_DUSP1_cyto'});
            testCase.Bursty.computeFIM([],'log');
        end
    end

    methods (TestMethodSetup)
        %
    end

    methods (Test)         
        function CalibratePDO(testCase)

            % Calibrate a PDO from average fluorescence intensity data
            testCase.ModelPDO_Intensity = testCase.Bursty.calibratePDO( ...
                'test_data/pdo.csv', ...
                {'rna'}, {'RNA_DUSP1_nuc'}, {'Nuc_DUSP1_avg_int'}, ...
                'AffinePoiss', true, [1, 1230, 0.5]);
        end

        function PadConditionalPmfs(testCase)

            % Calibrate PDOs for single Z slice
            testCase.ModelPDO_single_Z = testCase.Bursty.calibratePDO('test_data/Zslice.csv',...
                    {'rna'},{'Spot_Count'},{'Z_7'},'AffinePoiss',true);

            % Compute FIM for single Z slice
            fimsPDO_Z = testCase.ModelPDO_single_Z.computeFIM([],'log');

            % check
            testCase.verifyNotEmpty(fimsPDO_Z);
        end
    end
end
