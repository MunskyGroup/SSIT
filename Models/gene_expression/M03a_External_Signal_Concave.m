function concave_amplify_signal(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.
        TMP(1,1:5) = {'1','-','x1(1)','kx*I01','n'};
        TMP(2,:) = {'2','x1(1)','-','gx*x1','n'};
        TMP(3,:) = {'3','-','x2(1)','ky*x1/(1+x1)','n'};
        TMP(4,:) = {'4','x2(1)','-','gy*x2','n'};
        app.ReactionsTabOutputs.parameters(:,1) = {'kx','O','gx','ky','gy'};
        app.ReactionsTabOutputs.presetParameters = {10, 1, 10, 15, 1};
        app.ReactionsTabOutputs.inputs(:,1) = {'I01'};
        app.ReactionsTabOutputs.presetInputs = {'1+cos(O*t)'};
        app.ReactionsTabOutputs.initialCondition = '[0,0,0]';
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'Two species system with concave reaction rates to amplify external signals';'';'kx=production rate of x';'gx=degradation rate of species x';'ky=production rate of y';'gy=degradation rate of species y';'O=frequency of input signal'};
        app.ModelReactionTable.Data = TMP;
end