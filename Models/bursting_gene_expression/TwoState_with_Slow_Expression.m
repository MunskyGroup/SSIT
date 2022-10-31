function TwoState_with_Slow_Expression(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','x1(1)','x2(1)','kon_*(x1)','n'};
        TMP(2,:) = {'2','x2(1)','x1(1)','koff_*(x2)','n'};
        TMP(3,:) = {'3','-','x3(1)','krsl_*(x1)','n'};
        TMP(4,:) = {'4','-','x3(1)','krfs_*(x2)','n'};
        TMP(5,:) = {'5','x3(1)','-','kg_*(x3)','n'};
        app.ReactionsTabOutputs.presetParameters = {};
        app.ReactionsTabOutputs.presetInputs = {};
        app.ReactionsTabOutputs.initialCondition = '[1,0,0]';
        app.ReactionsTabOutputs.citations = {};
        app.ReactionsTabOutputs.modelInfo = {};
        app.ModelReactionTable.Data = TMP;
end