function Gene_Toggle_Model(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','-','x1(1)','ket*(ka1+((kb1*(k1^3))/((k1^3)+(x2)^3)))','n'};
        TMP(2,:) = {'2','x1(1)','-','(kd1+((ks*kg)/(1+ks)))*(x1)','n'};
        TMP(3,:) = {'3','-','x2(1)','ket*(ka2+((kb2*(k2^3))/((k2^3)+(x1)^3)))','n'};
        TMP(4,:) = {'4','x2(1)','-','kd2*(x2)','n'};
        scale = 20;
        app.ReactionsTabOutputs.parameters(:,1) = {'k1',    'k2',     'ka1',     'ka2',   'kb1',      'kb2', 'kd1','kd2', 'ks',  'kg',   'ket'};
        app.ReactionsTabOutputs.presetParameters = {1.0*scale,1.0*scale,0.2*scale,0.2*scale,4.0*scale,4.0*scale,1.0,  1.0,  1.0,  1.0,  0.1};
        app.ReactionsTabOutputs.presetInputs = {};
        app.ReactionsTabOutputs.initialCondition = '[0,0,0]';
        app.ReactionsTabOutputs.citations = {'T.S. Gardner, C.R. Cantor, and J.J. Collins. Construction of a genetic toggle switch in Escherichia coli. Nature, 403(6767):339�342, 2000.'};
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'genetic toggle switch in E.coli with two repressors and two promoters where each promoter is inhibited by the repressor that is transcribed by the opposing promoter'};
        app.ModelReactionTable.Data = TMP;
        app.FspPrintTimesField.Value = '[0:0.1:10]';
        app.PrintTimesEditField.Value = '[0:0.1:10]';
end