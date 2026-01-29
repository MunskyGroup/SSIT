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
        
        function test_ssa_tab_options(tc)
            % Change to SSA Tab
            tc.choose(tc.GUI.TabGroup,'Stochastic Simulation');
            % Run SSA
            tc.press(tc.GUI.SsaRunButton);
        end

        function test_fsp_tab_options(tc)
            % Change to FSP Tab
            tc.choose(tc.GUI.TabGroup,'Finite State Projection');
            tc.GUI.FspMarginalTimeCreateMovieCheckBox.Value = true;
            tc.GUI.FspMeanVarShowVarianceCheckBox.Value = true;
            tc.GUI.FspMeanVarShowOdeCheckBox.Value = true;
            % tc.GUI.FspMeshCheckBox = true;
            
            % Run FSP with various buttons
            Buttons = {'FspRunButtom','FspAddConstraintButton','FspRunButtom',...
                'FspUpdatePlotButton','FspPlotMarginalsButton','FspDefaultButton',...
                'FspRunButtom','FspPlotMarginalsOverTimeButton','FspMeanVarPlotButton'...
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

        function test_load_complex_model(tc)
            tc.choose(tc.GUI.TabGroup,'Model Loading and Building');
            tc.GUI = loadModelBP(tc.GUI, [], 'tests/test_data/GRDusp1ModelTestLibrary.mat');

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

            % Test that a hybrid model with data can be loaded into GUI.
            tc.choose(tc.GUI.ChooseSSITModel,'ModelDUSP1_100nM');
            % Test that the model can be saved and plots of fits can be
            % generated.

            tc.choose(tc.GUI.TabGroup,'Data Loading and Fitting');

            % Run Data Loading and Fitting with various buttons
            Buttons = {'SolveandPlotButton'};
            for iB = 1:length(Buttons)
                tc.press(tc.GUI.(Buttons{iB}));
            end
        end

        function test_hybrid_and_pdo_options(tc)
            % Switch to Distortion and Hybrid Modeling tab
            tc.choose(tc.GUI.TabGroup,'Distortion and Hybrid Modeling');

            % Change Distortion to Binomial and reset all PDO parameters.         
            tc.choose(tc.GUI.DistortionTypeDropDown,'None'); drawnow
            tc.choose(tc.GUI.DistortionTypeDropDown,'Poisson'); drawnow
            tc.choose(tc.GUI.DistortionTypeDropDown,'Binomial'); drawnow
            tc.GUI.SSITModel.pdoOptions.PDO = [];
            tc.GUI.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS1 = 0;
            tc.GUI.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS2 = 0;

            % Select 'rna' to generate PDO and make plots.
            tc.choose(tc.GUI.SpeciesDropDown,'rna')

        end
        % function test_hybrid_and_pdo_fitting(tc)
        %     %% Change back to data loading and fitting.
        %     tc.choose(tc.GUI.TabGroup,'Data Loading and Fitting');
        % 
        %     % Change fitting options to only 10 iterations.
        %     tc.GUI.DataLoadingAndFittingTabOutputs.fitOptions.props.Display = 'iter';
        %     tc.GUI.DataLoadingAndFittingTabOutputs.fitOptions.props.MaxIter = 10;
        % 
        %     % Run Data Loading and Fitting with various buttons
        %     Buttons = {'SolveandPlotButton','FitModelButton','SolveandPlotButton'};
        %     for iB = 1:length(Buttons)
        %         tc.press(tc.GUI.(Buttons{iB}));
        %     end
        % end

        function test_fsp_tab_hybrid_an(tc)
            % Change to FSP Tab
            tc.choose(tc.GUI.TabGroup,'Finite State Projection');
            tc.GUI.FspMarginalTimeCreateMovieCheckBox.Value = true;
            tc.GUI.FspMeanVarShowVarianceCheckBox.Value = true;
            tc.GUI.FspMeanVarShowOdeCheckBox.Value = true;
            % tc.GUI.FspMeshCheckBox = true;

            % Run FSP with various buttons
            Buttons = {'FspRunButtom','FspAddConstraintButton','FspRunButtom',...
                'FspUpdatePlotButton','FspPlotMarginalsButton','FspDefaultButton',...
                'FspRunButtom','FspPlotMarginalsOverTimeButton','FspMeanVarPlotButton'...
                };
            for iB = 1:length(Buttons)
                tc.press(tc.GUI.(Buttons{iB}));
            end
        end
    end
end