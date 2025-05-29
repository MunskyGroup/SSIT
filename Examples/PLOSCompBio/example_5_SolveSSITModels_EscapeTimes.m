%% example_5_SolveSSITModels_EscapeTimes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%% Solve a first-passage time problem (escape times)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
%clear
%close all
addpath(genpath('../../'));

% example_1_CreateSSITModels

% Load our models from example_1_CreateSSITModels and inspect them:
Model.summarizeModel
STL1.summarizeModel

% Set the times at which distributions will be computed:
Model.tSpan = linspace(0,20,200);
STL1.tSpan = linspace(0,20,200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Solve escape times for the bursting gene example model 
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model:
    % Create a copy of the bursting gene model:
    Model_escape = Model;
    
    %% Specify a boundary for the escape calculation
    % Calculate the time until the mRNA concentration reaches 50
    Model_escape.fspOptions.escapeSinks.f = {'offGene';'onGene'};
    Model_escape.fspOptions.verbose = false;
    Model_escape.fspOptions.escapeSinks.b = [0.5;0.5];
    Model_escape = Model_escape.formPropensitiesGeneral('Model_escape');
    [fspSoln_escape,Model_escape.fspOptions.bounds] = Model_escape.solve;
    Model_escape.makePlot(fspSoln_escape,'escapeTimes',[],[],10)
    
    %% Specify a boundary for the escape calculation
    % Calculate the time until the mRNA concentration reaches 50
    Model_escape.fspOptions.escapeSinks.f = {'mRNA'};
    Model_escape.fspOptions.verbose = false;
    Model_escape.fspOptions.escapeSinks.b = 0.5;
    Model_escape = Model_escape.formPropensitiesGeneral('Model_escape');
    [fspSoln_escape,Model_escape.fspOptions.bounds] = Model_escape.solve;
    Model_escape.makePlot(fspSoln_escape,'escapeTimes',[],[],10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve escape times for the time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1:
% Create a copy of the time-varying STL1 yeast model:
STL1_escape = STL1;
STL1_escape = STL1_escape.formPropensitiesGeneral('STL1_escape');

% Solve for the escape time:
STL1_escape.fspOptions.escapeSinks.f = {'offGene';'onGene'};
STL1_escape.fspOptions.escapeSinks.b = [0.5;0.5];
[STL1_fspSoln_escape,STL1_escape.fspOptions.bounds] = STL1_escape.solve;
STL1_escape.makePlot(STL1_fspSoln_escape,'escapeTimes',[],[],10)
% Note that with the decaying transcription rate not all cells will
% reach the level of 50 proteins.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve escape times for more complex escape thresholds.
%  In this example we explore the escape time until the number of proteins
%  proteins is more than 1.25 times the current number of mRNA molecules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Model_escape_complex = Model_escape;
Model_escape_complex.fspOptions.escapeSinks.f = {'mRNA/offGene'};
Model_escape_complex.fspOptions.escapeSinks.b = 1.25;
[fspSoln3,Model_escape_complex.fspOptions.bounds] = ...
                                                Model_escape_complex.solve;
Model_escape_complex.makePlot(fspSoln3,'escapeTimes',[],[],10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(4): Solve for the time until/probability that one specific condition
%% out of multiple possible escape conditions is met.
% In this example, we assume that there are two potential avenues to
% escape, and we want to know when and with what probability will each
% decision be made.  We want to know when:
%   A) the number of RNA exceeds 33 
%   B) the number of proteins exceeds 15
%   C) the number of proteins and RNA combined exceeds 45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Model_escape_decision = Model_escape_complex;
Model_escape_decision.fspOptions.escapeSinks.f = {'onGene';'mRNA';...
                                                  'onGene+mRNA'};
Model_escape_decision.fspOptions.escapeSinks.b = [0.5;0.5;1];
[fspSoln_esc_dec,Model_escape_decision.fspOptions.bounds] = ...
                                               Model_escape_decision.solve;
Model_escape_decision.makePlot(fspSoln_esc_dec,'escapeTimes',[],[],12)
legend(Model_escape_decision.fspOptions.escapeSinks.f)


