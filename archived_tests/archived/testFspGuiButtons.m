classdef test_FSP_GUI_buttons < matlab.uitest.TestCase
   properties
       app
   end
   methods (TestMethodSetup)
       function launchapp(tc)
           tc.app = SSA_and_FSP_Dan_and_Jim;
           tc.addTeardown(@delete,tc.app);
       end
   end
    methods (Test)
        
        function test_AddReactionButtonValueChanged(tc)
            tc.press(tc.app.AddReactionButton)
        end
        function test_RunFSPButtonPushed(tc)
            tc.press(tc.app.FspRunButton);
        end
        function test_UpdatePlotsButton_2(tc)
            tc.press(tc.app.FspUpdatePlotButton);
        end
        function test_x1CheckBox_2(tc)
            tc.app.Fspx1CheckBox.Value=0;
            tc.choose(tc.app.Fspx1CheckBox,true);
            tc.app.Fspx1CheckBox
            tc.verifyEqual(tc.app.Fspx1CheckBox.Value,false)
        end
        function test_x2CheckBox_2(tc)
            tc.press(tc.app.Fspx2CheckBox);
        end
        function test_x3CheckBox_2(tc)
            tc.press(tc.app.Fspx3CheckBox);
        end
        function test_SsaRunButton(tc)
            tc.press(tc.app.SsaRunButton);
            tc.app.Rslts
        end

    end
end