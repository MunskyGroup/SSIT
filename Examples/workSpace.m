%% Dusp1Models
clear all
clc
addpath('../CommandLine');
%% Define a "Complex" SSIT Model
% Here we define a model where there is a time varying concentration of
% Dexamethosone (IDex(t)).  Dex drives the transport of the GR to the
% nucleus at the rate (kcn0+kcn1*Dex), where 'kcn0' is the basal transport
% rate, and 'dcn1*IDex(t)' is the Dex-dependent rate. When in the nucleus,
% GR can activate an inactive gene with rate 'kon', and that can deactivate
% with rate 'koff'. Active genes can transcribe at rate 'kr' and the rna
% decay or move out of the nucleus with rate 'gr'.
% Here is how we create the whole model:
Model0 = SSIT;
Model0.species = {'offGene';'onGene';'cytGR';'nucGR';'rna'};
Model0.initialCondition = [2;0;0;40;5];
Model0.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
    '(kcn0+kcn1)*cytGR';'knc*nucGR';...
    'kr*onGene';'gr*rna'};
% Model0.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
%     '(kcn0+kcn1*IDex)*cytGR';'knc*nucGR';...
%     'kr*onGene';'gr*rna'};
Model0.stoichiometry = [-1,1,0,0,0,0;...
                         1,-1,0,0,0,0;...
                         0,0,-1,1,0,0;...
                         0,0,1,-1,0,0;...
                         0,0,0,0,1,-1];
% Model0.inputExpressions = {'IDex','exp(-gDex*t)*(t>0)'};
Model0.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
    'kcn0',0.0;'kcn1',0.0;'gDex',0.1;'knc',0.1});
% Model0.fspOptions.verbose = true;
Model0.fspOptions.initApproxSS = false;
Model0.tSpan = linspace(0,200,21);
[fspSoln0, Model0.fspOptions.bounds] = Model0.solve;
tic
[fspSoln0, Model0.fspOptions.bounds] = Model0.solve(fspSoln0.stateSpace);
toc
Model0.makePlot(fspSoln0,'marginals',[],[],[11:15])
Model0.makePlot(fspSoln0,'meansAndDevs',[],[],[1:5])

%% Assume that the GR signal dynamics are continuous and deterministic
% Here we show how to implement a simplified version of the model, where
% the GR signaling dynamics are assumed to evolve according to a simple ODE
% with the same stoichiometry and rates of the above system.  To achieve
% this, we simply copy the model and assign that the 'nucGR' species is
% deterministic.  We then solve and compare to the original solution.
Model2 = Model0;
Model2.useHybrid = true;
Model2.hybridOptions.upstreamODEs = {'nucGR','cytGR'};
[fspSoln2, Model2.fspOptions.bounds] = Model2.solve;
tic
Model2.fspOptions.bounds = Model0.fspOptions.bounds([1,2,5,6,7,10]);
[fspSoln2, Model2.fspOptions.bounds] = Model2.solve;
toc
Model2.makePlot(fspSoln2,'marginals',[],[],[11,12,15])
% For these parameters the approximate solution is very accurate, but is
% about 2x faster.

%%
Model0.fspOptions.verbose = false;
Model0.fspOptions.bounds(6) = 200; 
Model0 = Model0.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'rna','RNA_nuc'});
Model0.fittingOptions.modelVarsToFit = 1:7;
% Next, we call a fitting routine to maximize the likelihood of the data
% given the model.  Once that is complete, we update the model parameters
% and call a function to generate a plot of the results.
fitOptions = optimset('Display','iter','MaxIter',1000);
fitOptions.suppressFSPExpansion = true; 
% Fitting can be much faster if we choose not to expand the FSP during each
% step, but this also introduces an approximation error.
fitParameters = Model0.maximizeLikelihood([],fitOptions);
Model0.parameters(Model0.fittingOptions.modelVarsToFit,2) = num2cell(fitParameters);
Model0.makeFitPlot;

% Model1 = SSIT;
% Model1.species = {'onGene';'rna'};
% Model1.initialCondition = [0;0];
% Model1.propensityFunctions = {'kon*IGR*(2-onGene)';'koff*onGene';'kr*onGene';'gr*rna'};
% Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
% Model1.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
% Model1.parameters = ({'koff',0.14;'kon',0.14;'kr',1;'gr',0.01;...
%     'a1',0.4;'r1',0.04;'r2',0.1});
% Model1.fspOptions.initApproxSS = true;


%%
Model0 = SSIT;
Model0.species = {'rna'};
Model0.initialCondition = [0];
Model0.propensityFunctions = {'k*I';'gr*rna';};
Model0.stoichiometry = [1,-1];
Model0.inputExpressions = {'I','(t<1)'};
Model0.tSpan = linspace(0,10,21);
Model0.parameters = ({'k',50;'gr',1});
Model0.fspOptions.initApproxSS = false;
tic
[fspSoln0, Model0.fspOptions.bounds] = Model0.solve;
toc
Model0.makePlot(fspSoln0,'meansAndDevs',[],[])

