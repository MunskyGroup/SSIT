classdef testGUI < matlab.uitest.TestCase
    properties
        GUI
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function createTestGUI(tc)
            cd ../../SSIT/
            addpath('../')
            addpath(genpath('../src'));
            close all force
            tc.GUI = SSIT_GUI;
        end
    end

    methods (TestMethodSetup)

    end

    methods (Test)
        % Test methods

        function GUIExists(tc)
            tc.verifyEqual(isa(tc.GUI,'SSIT_GUI'), true, ...
                'the SSIT GUI was launched');
        end

        function ModelUpdateButton(tc)
            tc.choose(tc.GUI.ModelDropDown,'M01_Toggle_Switch.m')
            tc.verifyEqual(strcmp(tc.GUI.ModelDropDown.Value,...
                'M01_Toggle_Switch.m'), true, ...
                'The SSIT model was changed');
        end

        function runSSA(tc)
            tc.choose(tc.GUI.StochasticSimulationTab)
            tc.press(tc.GUI.SsaRunButton)
        end

        function runFSP(tc)
            tc.choose(tc.GUI.FSPTab)
            tc.press(tc.GUI.FspRunButtom)
        end

        function runSens(tc)
            tc.choose(tc.GUI.SensitivityTab)
            tc.press(tc.GUI.SensRunButton)
        end

        function runFIM(tc)
            tc.choose(tc.GUI.FisherInformationTab)
            tc.press(tc.GUI.PlotInformationvsTimeButton)
        end

    end

end