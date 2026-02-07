function [] = createNewConstraint(app)
% Adds a new constraint line to the Constraints table within the FSP tab of
% the GUI
TMP = app.FspConstraintTable.Data;

speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);
nSpecies = length(speciesStochastic);
           
tmp = ['(',speciesStochastic{1},'-3)'];
for i = 2:nSpecies
    tmp = [tmp,'*(',speciesStochastic{i},'-3)'];
end
TMP{end+1,1} = tmp;
TMP{end,2} = '<';
TMP{end,3} = 1;
app.FspConstraintTable.Data=TMP;
app.FspTabOutputs.stateSpace = [];
if length(app.SSITModel.fspOptions.bounds)<size(TMP,1)
    app.SSITModel.fspOptions.bounds(length(app.SSITModel.fspOptions.bounds)+1:size(TMP,1)) = ...
        [TMP{length(app.SSITModel.fspOptions.bounds)+1:size(TMP,1),3}];
end
end