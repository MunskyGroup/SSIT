function toggle_switch(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','-','x1(1)','kn0+kn1*(1/(1+asn*(x2^nsn)))','n'};
        TMP(2,1:5) = {'2','-','x2(1)','kn0+kn1*(1/(1+asn*(x1^nsn)))','n'};
        TMP(3,1:5) = {'3','x1(1)','-','g*x1','n'};
        TMP(4,1:5) = {'4','x2(1)','-','g*x2','n'};
        app.ReactionsTabOutputs.parameters(:,1) = {'kn0','kn1','asn','nsn','g'};
        app.ReactionsTabOutputs.presetParameters = {1, 50, 5, 2, 1};
        app.ReactionsTabOutputs.presetInputs = {};
        app.ReactionsTabOutputs.initialCondition = '[0,0,0]';
        app.ReactionsTabOutputs.citations = {'Tapia, J., Faeder, J., Munsky, B., Adaptive Coarse-Graining for Transient and Quasi-Equilibrium Analyses of Stochastic Gene Regulation, Proc. of the 51st IEEE Conference on Decision and Control, 5361-5366, Maui, HI, Dec. 2012.'};
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'Genetic toggle switch in which two genes are negatively regulated by each other';'';'asn=ratio of kon to koff';'kn0=basal production rate';'kn1=active production rate';'nsn=cooperativity of repression';'g=degradation rate'};
        app.ModelReactionTable.Data = TMP;
end