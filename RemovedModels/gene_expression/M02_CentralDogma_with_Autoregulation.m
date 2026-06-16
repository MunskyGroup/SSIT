function autoregulation_feedback(app)
% This function defines the model parameters for the FSP 3 Species GUI.
% It fills in the reactions table in the reactions tab.
        TMP(1,1:5) = {'1','-','x1(1)','kr/(1+kf*x2)','n'};
        TMP(2,:) = {'2','x1(1)','-','gr*x1','n'};
        TMP(3,:) = {'3','-','x2(1)','kp*x1','n'};
        TMP(4,:) = {'4','x2(1)','-','gp*x2','n'};
        app.ReactionsTabOutputs.parameters(:,1) = {'kr','kf','gr','kp','gp'};
        app.ReactionsTabOutputs.presetParameters = {120, 1, 1, 5, 1};
        app.ReactionsTabOutputs.presetInputs = {};
        app.ReactionsTabOutputs.initialCondition = '[0,0,0]';
        app.ReactionsTabOutputs.citations = {'Munsky, B., Modeling Cellular Variability, in Quantitative Biology From Molecular to Cellular Systems, pp. 234-266, M. Wall, Ed. (Taylor & Francis Group, New York, 2012).https://www.engr.colostate.edu/~munsky/Papers/Munsky_Chapter_2012.pdf'};
        app.ReactionsTabOutputs.modelInfo={'About the Model';'';'Autoregulation of a single gene whose protein inhibits its own transcription';'';'kf=strength of negative feedback';'kr=basal transcription rate';'gr=degradation rate of mRNA transcripts';'kp=production rate of the protein';'gp= degradation rate of proten'};
        app.ModelReactionTable.Data = TMP;
end
