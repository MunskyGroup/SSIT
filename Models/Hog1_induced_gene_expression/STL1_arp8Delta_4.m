function STL1_arp8Delta_4(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','-','x1(1)','k12*(x1==1)+k23*(x1==2)+k34*(x1==3)','n'};
        TMP(2,:) = {'2','x1(1)','-','k21*(x1==2)*IHog+k32*(x1==3)+k43*(x1==4)','n'};
        TMP(3,:) = {'3','-','x2(1)','kr1*(x1==1)+kr2*(x1==2)+kr3*(x1==3)+kr4*(x1==4)','n'};
        TMP(4,:) = {'4','x2(1)','-','kg*(x2)','n'};
        app.ReactionsTabOutputs.parameters(:,1) = {'k12','k23','k34','k21','k32','k43','kr1','kr2','kr3','kr4','kg'};
        app.ReactionsTabOutputs.presetParameters = {0.0039,0.0013, 0.16, 1, 0.021, 0.038, 1.7e-4, 0.0020, 1, 0.054, 0.0049};
        app.ReactionsTabOutputs.inputs(:,1) = {'IHog'};
        app.ReactionsTabOutputs.presetInputs = {'max(0,120*(1-(2.4*(9.3e9*(((1-exp(-(2.2e-9)*t))*exp(-(2.6e-3)*t))/(1+(((1-exp(-(2.2e-9)*t))*exp(-(2.6e-3)*t))/(4.1e-7))))^1.6))))'};
        app.ReactionsTabOutputs.initialCondition = '[1,0,0]';
        app.ReactionsTabOutputs.citations = {'Neuert, G., Munsky, B., Tan, R. Z., Teytelman, L., Khammash, M., & van Oudenaarden, A. (2013). Systematic Identification of Signal-Activated Stochastic Gene Regulation. Science, 339(6119), 584�587.'};
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'high-osmolarity glycerol (HOG) mitogen-activated protein kinase (MAPK) pathway in S. cerevisiae and its regulation of the STL1 arp8Delta mutant gene through osmotic shock.';'';'kij=conversion rate from state i to j';'kri=production rate of mRNA from state i';'kg=degradation rate'};
        app.ModelReactionTable.Data = TMP;
        
end