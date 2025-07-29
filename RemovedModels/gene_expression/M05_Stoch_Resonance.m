function resonance(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.
        TMP(1,1:5) = {'1','-','x1(1)','((x1==1)*x2^10)/(k1^10+x2^10)','n'};
        TMP(2,:) = {'2','x1(1)','-','(k2*(x1==2))/(1+x2^10)','n'};
        TMP(3,:) = {'3','-','x2(1)','k3*(x1==1)','n'};
        TMP(4,:) = {'4','x2(1)','-','k4*x2','n'};
        app.ReactionsTabOutputs.parameters(:,1) = {'k1','k2','k3','k4'};
        app.ReactionsTabOutputs.presetParameters = {1000, 100, 20, 0.02};
        app.ReactionsTabOutputs.presetInputs = {};        
        app.ReactionsTabOutputs.citations = {'Munsky, B., Modeling Cellular Variability, in Quantitative Biology From Molecular to Cellular Systems, pp. 234-266, M. Wall, Ed. (Taylor & Francis Group, New York, 2012).https://www.engr.colostate.edu/~munsky/Papers/Munsky_Chapter_2012.pdf'};
        app.ReactionsTabOutputs.initialCondition = '[1,100,0]';
        app.ModelReactionTable.Data = TMP;
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'Theoretical model of circadian rhythm with a single gene with two state (s1 or s2) where Y is produced in the s1 state';'';'k1=state transition rate to s2';'k2=state transition rate to s1';'k3=degradation rate of s2';'k4=degradation of y'};
end