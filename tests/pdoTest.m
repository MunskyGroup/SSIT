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

            % Calibrate PDOs for single and subset Z slices
            testCase.ModelPDO_single_Z = testCase.Bursty.calibratePDO('../WorkSpace/AlexP/SpotCountingByCellAndByZ/merged_dataframe_A549_DUSP1_100nM_10min_062723.csv',...
                    {'rna'},{'Spot_Count'},{'Z_7'},'AffinePoiss',true);
            testCase.ModelPDO_5thru9_Zs = testCase.Bursty.calibratePDO('../WorkSpace/AlexP/SpotCountingByCellAndByZ/merged_dataframe_A549_DUSP1_100nM_10min_062723.csv',...
                    {'rna'},{'Spot_Count'},{'Z_5_9'},'AffinePoiss',true);

            % Compute FIM for single Z and Zs 5-9
            fimsPDO_Z = testCase.ModelPDO_single_Z.computeFIM([],'log');
            fimsPDO_5thru9_Zs = testCase.ModelPDO_5thru9_Zs.computeFIM([],'log');

            % check
            testCase.verifyNotEmpty(fimsPDO_Z);
            testCase.verifyNotEmpty(fimsPDO_5thru9_Zs);
        end
    end
end
