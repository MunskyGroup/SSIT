clear all
clc

%% Complex Dusp1 model
Model = SSIT;
Model.species = {'x1';'x2';'x3'};  % GRnuc, geneOn, dusp1
Model.initialCondition = [0;0;0];
Model.propensityFunctions = {'(kcn0+kcn1*IDex)';'knc*x1';...
  'kon*x1*(2-x2)';'koff*x2';'kr*x2';'gr*x3'};
Model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'};
Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
  'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
Model.stoichiometry = [ 1,-1, 0, 0, 0, 0;...
                      0, 0, 1,-1, 0, 0;...
                      0, 0, 0, 0, 1,-1];
Model.fspOptions.initApproxSS = true;
%% load data
mhResults = load('complex_dusp1_mhast.mat').mhResults;
sensSoln = load('complex_dusp1_sens.mat').sensSoln;
comp_Model = load('complex_dusp1_model.mat').Model;
%% Tryptolide Experiment
ModelTrypt = comp_Model;
ModelTrypt.propensityFunctions(5) = {'kr*x2*Itrypt'};
ModelTrypt.inputExpressions(2,:) = {'Itrypt','(t<tpt)'};


ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';

tpt_array = 55;
% sensSoln = cell(1,length(tpt_array));
for iTpt = 1:length(tpt_array)
    iTpt
    ModelTrypt.parameters(9,:) = {'tpt',tpt_array(iTpt)};

    % Get FSP fit for bounds.
    ModelTrypt.solutionScheme = 'FSP';
    ModelTrypt.fspOptions.fspTol = 1e-6;
    ModelTrypt.fspOptions.bounds=[];
    [fspSolntpt,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;

    ModelTrypt.fspOptions.fspTol = inf;
    ModelTrypt.solutionScheme = 'fspSens';
%     sensSoln{iTpt} = ModelTrypt.solve(fspSoln.stateSpace);

    % Solve for Fisher Information Matrix at all Time Points
%     ModelTrypt.pdoOptions.unobservedSpecies = {'x1','x2'};
%     fims = ModelTrypt.computeFIM(sensSoln{iTpt}.sens);
%     FIM = ModelTrypt.evaluateExperiment(fims,ModelTrypt.dataSet.nCells);
%     expectedDetCov(iTpt) = det(FIM^(-1));
end
%%
figure;plot(tpt_array,expectedDetCov); hold on
set(gca,'yscale','log')
ylim = get(gca,'ylim');
for it = 1:length(ModelTrypt.dataSet.times)
    plot(ModelTrypt.dataSet.times(it)*[1,1],ylim,'k--')
end