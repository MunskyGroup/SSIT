function STL1_WT_6(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.

        TMP(1,1:5) = {'1','-','x1(1)','k12*(x1==1)+k23*(x1==2)+k34*(x1==3)','n'};
        TMP(2,:) = {'2','x1(1)','-','k21*(x1==2)*I1+k32*(x1==3)+k43*(x1==4)','n'};
        TMP(3,:) = {'3','-','x2(1)','kr1*(x1==1)+kr2*(x1==2)+kr3*(x1==3)+kr4*(x1==4)','n'};
        TMP(4,:) = {'4','x2(1)','-','g*(x2)','n'};
        app.ReactionsTabOutputs.parameters(:,1) = {'k12','k21','k23','k32','k34','k43','kr1','kr2','kr3','kr4','g'};
        app.ReactionsTabOutputs.presetParameters = {1.3,1,0.0067, 0.027, 0.13, 0.038, 7.8e-4, 0.012, 0.99, 0.0540, 0.0049};
        app.ReactionsTabOutputs.inputs(:,1) = {'I1'};
        app.ReactionsTabOutputs.presetInputs = {'max(0,3200*(1-(2.4*(9.3e9*(((1-exp(-(6.9e-5)*t))*exp(-(2.4e-3)*t))/(1+(((1-exp(-(6.9e-5)*t))*exp(-(2.4e-3)*t))/(6.4e-4))))^3.1))))'};
        app.ReactionsTabOutputs.initialCondition = '[1,0,0]';
        app.ReactionsTabOutputs.citations = {'Neuert, G., Munsky, B., Tan, R. Z., Teytelman, L., Khammash, M., & van Oudenaarden, A. (2013). Systematic Identification of Signal-Activated Stochastic Gene Regulation. Science, 339(6119), 584�587.'};
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'high-osmolarity glycerol (HOG) mitogen-activated protein kinase (MAPK) pathway in S. cerevisiae and its regulation of the STL1 gene through different levels of osmotic shock.';'';'kij=conversion rate from state i to j';'kri=production rate of mRNA from state i';'kg=degradation rate'};
        app.ModelReactionTable.Data = TMP;
        
end