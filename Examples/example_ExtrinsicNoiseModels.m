%% example_ExtrinsicNoiseModels
% In this example, we show how to sample an FSM  model over intrinsic noise
% in its various parameters. 
close all 
clear all
addpath(genpath('../src'));
%% Example 1 - transcription and translation
% First create a full model (e.g., for mRNA and protein)
Model1 = SSIT();
Model1.species = {'rna','protein'};
Model1.initialCondition = [0;0];
Model1.propensityFunctions = {'kr';'gr*rna';'k2*rna';'g2*protein'};
Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model1.parameters = ({'kr',10;'gr',0.5;...
    'k2',2;'g2',1});
Model1.fspOptions.initApproxSS = false; 
Model1.tSpan = linspace(0,5,10);
[fspSoln1,Model1.fspOptions.bounds] = Model1.solve;
Model1 = Model1.formPropensitiesGeneral('Model1');   
Model1.makePlot(fspSoln1,'marginals',[],[],[2,3])

%% Generate, solve and plot results for an extrinsic noise version
% Specify the rules for the extrinsic noise.  This must be a function that
% returns a  vector of parameters, that are in the same order as the
% parameters provided. You can choose a different distribution for the
% extrinsic noise in each parameter.
parDistributions = @()[10+2*randn,...
    0.5+0.1*randn,...
    2+0.4*randn,...
    1+0.2*randn];
extrinsicModel = extrinsicSSIT(Model1,parDistributions,20);
Model1.makePlot(extrinsicModel.averagedResults,'marginals',[],[],[2,3])