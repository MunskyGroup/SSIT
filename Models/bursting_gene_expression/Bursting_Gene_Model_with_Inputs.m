function Bursting_Gene_Model_with_Inputs(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','-','x1(1)','kon_*(x1==0)*I1_','n'};
        TMP(2,:) = {'2','x1(1)','-','koff_*(x1==1)','n'};
        TMP(3,:) = {'3','-','x2(1)','kr_*(x1==1)','n'};
        TMP(4,:) = {'4','x2(1)','-','kg_*(x2)','n'};
        TMP(5,:) = {'5','-','x3(1)','kp_*(x2)','n'};
        TMP(6,:) = {'6','x3(1)','-','kgp_*(x3)','n'};
        app.ReactionsTabOutputs.presetParameters = {};
        app.ReactionsTabOutputs.presetInputs = {};
        app.ReactionsTabOutputs.initialCondition = '[0,0,0]';
        app.ReactionsTabOutputs.citations = {};
        app.ReactionsTabOutputs.modelInfo = {};
        app.ModelReactionTable.Data = TMP;
end