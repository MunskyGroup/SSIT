classdef testGui < matlab.uitest.TestCase
   properties
       GUI
   end
   methods (TestClassSetup)
       function launchapp(testCase1)

           addpath(genpath('../src'));
           
           stubDir = tempname;
           mkdir(stubDir);

           % --- Write stub questdlg to that folder ---
           fid = fopen(fullfile(stubDir, 'questdlg.m'), 'w');
           fprintf(fid, 'function sel = questdlg(varargin)\n');
           fprintf(fid, 'sel = ''Yes - Overwrite'';\n');  % always return "Yes"
           fprintf(fid, 'end\n');
           fclose(fid);

           % --- Put stub in front of MATLAB path ---
           origPath = addpath(stubDir);
           cleanup = onCleanup(@() path(origPath)); % restore later

           close all force
           cd ../  % Change to the main SSIT directory.
           testCase1.GUI = SSITGUI;
           testCase1.addTeardown(@delete,testCase1.GUI);
       end
   end
    methods (Test)
        function test_change_model_folder(tc)
            tc.press(tc.GUI.RefreshFoldersButton)
        end            
        
        function test_select_model(tc)
            tc.choose(tc.GUI.ModelDropDown,'M00_Poisson_Process.mat');
        end
        
        function test_change_to_ssa_tab(tc)
            % Change to SSA Tab
            tc.choose(tc.GUI.TabGroup,'Stochastic Simulation');
            % Run SSA
            tc.press(tc.GUI.SsaRunButton);
        end

        function test_change_to_fsp_tab(tc)
            % Change to FSP Tab
            tc.choose(tc.GUI.TabGroup,'Finite State Projection');
            tc.GUI.FspMarginalTimeCreateMovieCheckBox.Value = true;
            tc.GUI.FspMeanVarShowVarianceCheckBox.Value = true;
            tc.GUI.FspMeanVarShowOdeCheckBox.Value = true;
            % tc.GUI.FspMeshCheckBox = true;
            
            % Run FSP with various buttons
            Buttons = {'FspRunButtom','FspAddConstraintButton','FspRunButtom',...
                'FspUpdatePlotButton','FspPlotMarginalsButton','FspDefaultButton',...
                'FspRunButtom','FspPlotMarginalsOverTimeButton','FspMeanVarPlotButton',...
                };
            for iB = 1:length(Buttons)
                tc.press(tc.GUI.(Buttons{iB}));
            end
        end
        
        function test_change_to_sens_tab(tc)
            % Change to Sensitivity Tab
            tc.choose(tc.GUI.TabGroup,'Sensitivity Analysis');

            % Run FSP with various buttons
            Buttons = {'SensRunButton'};
            for iB = 1:length(Buttons)
                tc.press(tc.GUI.(Buttons{iB}));
            end
        end

        function test_change_to_FIM_tab(tc)
            % Change to Sensitivity Tab
            tc.choose(tc.GUI.TabGroup,'Fisher Information');

            hold(tc.GUI.FIMEllipseAxes,'on');

            % Run FSP with various buttons
            Buttons = {'PlotInformationvsTimeButton',...
                'ManuallyAllocateMeasurementsperTimePointButton',...
                'EstimateMLEUncertaintyButton',...
                'OptimizeButton'};
            for iB = 1:length(Buttons)
                tc.press(tc.GUI.(Buttons{iB}));
            end
        end
    end
end