%% example_3_EscapeTimes
% Example script to demonstrate how to create and solve a first-passage 
% time problem.
close all
clear 
addpath(genpath('../../src'));

%% Preliminaries
% Load our models from example_1_CreateSSITModels and inspect them:
example_1_CreateSSITModels
Model.summarizeModel
STL1Model.summarizeModel

% Set the times at distributions will be computed:
Model.tSpan = linspace(0,20,200);
STL1Model.tSpan = linspace(0,20,200);

Model_escape = Model;
STL1Model_escape = STL1Model;

%% Ex.(1) 
%% Example 1 - a simple transcription/translation model
% First create a full model (e.g., for mRNA and protein)
% Model1 = SSIT;
% Model1.species = {'rna','protein'};
% Model1.initialCondition = [0;0];
% Model1.propensityFunctions = {'kr';'gr*rna';'k2*rna';'g2*protein'};
% Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
% Model1.parameters = ({'kr',100;'gr',0.5;...
%     'k2',2;'g2',1});

%% Specify a boundary for the escape calculation
% Calculate the time until the mRNA concentration reaches 50.
Model_escape.fspOptions.escapeSinks.f = {'offGene';'onGene'};
Model_escape.fspOptions.verbose = false;
Model_escape.fspOptions.escapeSinks.b = [0.5;0.5];
Model_escape = Model_escape.formPropensitiesGeneral('Model_escape');
[fspSoln_escape,Model_escape.fspOptions.bounds] = Model_escape.solve;
Model_escape.makePlot(fspSoln_escape,'escapeTimes',[],[],10)

%% Specify a boundary for the escape calculation
% Calculate the time until the mRNA concentration reaches 50.
Model_escape.fspOptions.escapeSinks.f = {'mRNA'};
Model_escape.fspOptions.verbose = false;
Model_escape.fspOptions.escapeSinks.b = 0.5;
Model_escape = Model_escape.formPropensitiesGeneral('Model_escape');
[fspSoln_escape,Model_escape.fspOptions.bounds] = Model_escape.solve;
Model_escape.makePlot(fspSoln_escape,'escapeTimes',[],[],10)

%% Example 2 - escape time with time varying transcription rate

STL1Model_escape = STL1Model_escape.formPropensitiesGeneral('STL1Model_escape');

% Solve for the escape time:
STL1Model_escape.fspOptions.escapeSinks.f = {'offGene';'onGene'};
STL1Model_escape.fspOptions.escapeSinks.b = [0.5;0.5];
[STL1Model_fspSoln_escape,STL1Model_escape.fspOptions.bounds] = STL1Model_escape.solve;
STL1Model_escape.makePlot(STL1Model_fspSoln_escape,'escapeTimes',[],[],10)
% Note that with the decaying transcription rate not all cells will
% reach the level of 50 proteins.

%% Example 3 - More complex escape thresholds.
% In this example we explore the escape time until the number of proteins
% is more than 1.25 times the current number of mRNA molecules.
Model_escape_complex = Model_escape;
Model_escape_complex.fspOptions.escapeSinks.f = {'mRNA/offGene'};
Model_escape_complex.fspOptions.escapeSinks.b = 1.25;
[fspSoln3,Model_escape_complex.fspOptions.bounds] = Model_escape_complex.solve;
Model_escape_complex.makePlot(fspSoln3,'escapeTimes',[],[],10)

%% Example 4 - Time until making a specific decision.
% In this example, we assume that there are two potential avenues to
% escape, and we want to know when and with what probability will each
% decision be made.
% In our case we want to know when:
%   A) the number of RNA exceeds 33 
%   B) the number of proteins exceeds 15
%   C) the number of proteins and RNA combined exceeds 45
Model_escape_decision = Model_escape_complex;
Model_escape_decision.fspOptions.escapeSinks.f = {'onGene';'mRNA';'onGene+mRNA'};
Model_escape_decision.fspOptions.escapeSinks.b = [0.5;0.5;1];
[fspSoln_esc_dec,Model_escape_decision.fspOptions.bounds] = Model_escape_decision.solve;
Model_escape_decision.makePlot(fspSoln_esc_dec,'escapeTimes',[],[],12)
legend(Model_escape_decision.fspOptions.escapeSinks.f)


