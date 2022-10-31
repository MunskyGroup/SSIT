function Five_Gene_State_Web_Model_with_2_RNAs_Expressed(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','-','x1(4)','k15_*(x1==1)','n'};
        TMP(2,:) = {'2','x1(3)','-','k52_*(x1==5)','n'};
        TMP(3,:) = {'3','-','x1(3)','k14_*(x1==1)','n'};
        TMP(4,:) = {'4','x1(1)','-','k21_*(x1==2) + k32_*(x1==3) + k43_*(x1==4)','n'};
        TMP(5,:) = {'5','x1(2)','-','k31_*(x1==3)','n'};
        TMP(6,:) = {'6','-','x2(1)','k2x2_*(x1==2)','n'};
        TMP(7,:) = {'7','-','x3(1)','k4x3_*(x1==4)','n'};
        TMP(8,:) = {'8','x1(3)','x3(1)','k5x3_*(x1==5)','n'};
        app.ReactionsTabOutputs.presetParameters = {};
        app.ReactionsTabOutputs.presetInputs = {};
        app.ReactionsTabOutputs.initialCondition = '[1,0,0]';
        app.ReactionsTabOutputs.citations = {};
        app.ReactionsTabOutputs.modelInfo = {};
        app.ModelReactionTable.Data = TMP;
end