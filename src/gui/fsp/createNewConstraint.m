function [] = createNewConstraint(app)
% Adds a new constraint line to the Constraints table within the FSP tab of
% the GUI
TMP = app.FspConstraintTable.Data;
nSpecies = length(app.SSITModel.species);
           
tmp = ['(',app.SSITModel.species{1},'-3)'];
for i = 2:nSpecies
    tmp = [tmp,'*(',app.SSITModel.species{i},'-3)'];
end
TMP{end+1,1} = tmp;
TMP{end,2} = '<';
TMP{end,3} = 1;
app.FspConstraintTable.Data=TMP;
app.FspTabOutputs.stateSpace = [];
end