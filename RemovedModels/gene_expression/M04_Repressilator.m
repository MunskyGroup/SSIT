function reaction_channels_repressilator(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','-','x1(1)','kn0+kn1*(1/(1+a*(x2^n)))','n'};
        TMP(2,1:5) = {'2','-','x2(1)','kn0+kn1*(1/(1+a*(x3^n)))','n'};
        TMP(3,1:5) = {'3','-','x3(1)','kn0+kn1*(1/(1+a*(x1^n)))','n'};
        TMP(4,1:5) = {'4','x1(1)','-','g*x1','n'};
        TMP(5,1:5) = {'5','x2(1)','-','g*x2','n'};
        TMP(6,1:5) = {'6','x3(1)','-','g*x3','n'};
        app.ReactionsTabOutputs.parameters(:,1) = {'kn0','kn1','a','n','g'};
        app.ReactionsTabOutputs.presetParameters = {0, 25, 5, 6, 1};
        app.ReactionsTabOutputs.presetInputs = {};
        app.ReactionsTabOutputs.initialCondition = '[20,0,0]';
        app.ReactionsTabOutputs.citations = {'Tapia, J., Faeder, J., Munsky, B., Adaptive Coarse-Graining for Transient and Quasi-Equilibrium Analyses of Stochastic Gene Regulation, Proc. of the 51st IEEE Conference on Decision and Control, 5361-5366, Maui, HI, Dec. 2012.'};
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'System of three chemical species that regulate each other through repression in a sequential feedback loop';'';'a=ratio of kon to koff';'kn0=basal production rate';'kn1=active production rate';'n=cooperativity of repression';'g=degradation rate'};
        app.ModelReactionTable.Data = TMP;
end