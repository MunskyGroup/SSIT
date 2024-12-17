classdef pdoTest < burstyTest
    properties
        
    end

    methods (TestClassSetup)
        % Define SSIT model	
         function setupModel(testCase)
             [testCase.Bursty, testCase.BurstyFSPSoln] = modelBuilder.buildBurstyModel();
         end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)  
        
        function CalibratePDO(testCase)
        	% Check PDO calibrations.
        	% Calibrate a PDO for Single Z slice data.
        	ModelPDO_Single_Z = testCase.Bursty.calibratePDO('../tests/test_data/pdo.csv',...
    			{'rna'},{'Spot_Count'},{'Z_7'},'AffinePoiss',true);
			% Calibrate a PDO from average fluorescence intensity data.
			ModelPDO_Intensity = testCase.Bursty.calibratePDO('../tests/test_data/pdo.csv',...
    			{'rna'},{'RNA_DUSP1_cyto'},{'Cyto_DUSP1_avg_int'},'AffinePoiss',true,[1,230,0.5]);
    	end
        
        % function ComputingLikelihood(testCase)               
        % end

        % function ComputingSensitivities(testCase)         
        % end
    % 
    %     function TestFIM(testCase)                 
    %     end
    % 
    %     function TestFIMwPrior(testCase)          
    %     end
    % 
    %     function FitUsingFSP(testCase)            
    %     end  
    % 
    %     function MetHastAndSampledFIM(testCase)
    %     end  
    % end
    end
end