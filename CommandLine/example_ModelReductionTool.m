clear all
clc

%% Define SSIT Model
% SSIT models are defined as usual:
% Model1 = SSIT;
% Model1.species = {'x1';'x2'};
% Model1.initialCondition = [0;0];
% Model1.propensityFunctions = {'kon*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
% Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
% Model1.parameters = ({'koff',14;'kon',14;'kr',25;'gr',0.1});
% Model1.fspOptions.initApproxSS = false;
% Model1.tSpan = linspace(0,10,1000);
% 
Model1 = SSIT;
Model1.species = {'x1';'x2'};
Model1.initialCondition = [0;0];
Model1.propensityFunctions = {'kr';'gr*x1';'kp*x1';'gp*x2'};
Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model1.parameters = ({'kr',20;'gr',1;'kp',1;'gp',1});
Model1.fspOptions.initApproxSS = false;
Model1.tSpan = linspace(0,5,1000);
% 
% Model1 = SSIT;
% Model1.species = {'x1'};
% Model1.initialCondition = [0];
% Model1.propensityFunctions = {'kr';'gr*x1'};
% Model1.stoichiometry = [1,-1];
% Model1.parameters = ({'kr',200;'gr',1});
% Model1.fspOptions.initApproxSS = false;
% Model1.tSpan = linspace(0,5,1000);

%%
% Model1.modelReductionOptions.useModReduction = false;
[fspSoln,Model1.fspOptions.bounds] = Model1.solve;
Model1.makePlot(fspSoln,'meansAndDevs',[],[],1)
Model1.makePlot(fspSoln,'marginals',[1,10,300,1000],[],[2,3])
%%
Model2 = Model1;
Model2.fspOptions.fspTol = inf;       
% Model2.modelReductionOptions.reductionType = 'Proper Orthogonal Decomposition';
% Model2.modelReductionOptions.reductionType = 'Dynamic Mode Decomposition';
% Model2.modelReductionOptions.reductionType = 'Eigen Decomposition Initial';
Model2.modelReductionOptions.reductionType = 'Logarithmic State Lumping';
Model2.modelReductionOptions.reductionOrder = 30;
Model2 = Model2.computeModelReductionTransformMatrices(fspSoln);

%%
Model2.modelReductionOptions.useModReduction = true;
fspSoln2 = Model2.solve(fspSoln.stateSpace);
Model2.makePlot(fspSoln2,'meansAndDevs',[],[],1)
Model2.makePlot(fspSoln2,'marginals',[1,10,300,1000],[],[2,3])
