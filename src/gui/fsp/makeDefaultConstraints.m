function [constr] = makeDefaultConstraints(app)
% This fills in the default constraints table
arguments
    app = [];
end

% if isempty(app)
%     app.ReactionsTabOutputs.varNames = obj.species;    
% end

Data={};

speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);
nSpecies = length(speciesStochastic);

for i = 1:nSpecies
    Data(i,1) = {['-',speciesStochastic{i}]};
    Data(i,2:3) = {'<',0};
    Data(nSpecies+i,1) = {speciesStochastic{i}};
    Data(nSpecies+i,2:3) = {'<',1};
end
app.FspConstraintTable.Data=Data;

% if nargout>=1
%     constr.f = Data(:,1);
%     constr.b = [Data{:,3}];
% end

app.FspTabOutputs.stateSpace = [];