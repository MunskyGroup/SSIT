%% example_2_SolveSSITModels_Hybrid
% Example script to demonstrate how to adjust a model to treat some species
% species (i.e., upstream reactions) using an ODE formulation, while having
% other species (i.e., downstream species) evolving in a discrete
% stochastic manner. This runs signficantly faster than full FSP solutions.
close all 
clear 
addpath(genpath('../../'));

%% Preliminaries
% Load our models from example_1_CreateSSITModels and inspect them:
example_1_CreateSSITModels
Model.summarizeModel

% Set the times at distributions will be computed:
Model_hybrid.tSpan = linspace(0,20,200);

%% Compute Ordinary Differential Equations (ODEs)

    % Create a copy of Model
    Model_hybrid = Model;

    % Set 'useHybrid' to true
    Model_hybrid.useHybrid = true;

    % Define which species will be solved by ODEs
    Model_hybrid.hybridOptions.upstreamODEs = {'offGene','onGene'};

    % Set solution scheme to FSP
    Model_hybrid.solutionScheme = 'FSP';

    % View summary of hybrid model, expecting:
    %   Species:
    %       offGene; IC = 1;  upstream ODE
    %       onGene; IC = 0;  upstream ODE
    %       mRNA; IC = 0;  discrete stochastic
    Model_hybrid.summarizeModel

    % Optionally, define a custom constraint on the species 
    % (e.g.,'offGene'+'onGene')
    Model_hybrid.customConstraintFuns = [];

    % Set FSP 1-norm error tolerance:
    Model_hybrid.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    Model_hybrid.fspOptions.bounds = [2,2,400];

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'Model_ODE'
    Model_hybrid = Model_hybrid.formPropensitiesGeneral('Model_hybrid');

    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    Model_hybrid.fspOptions.initApproxSS = false; 
    
    % Solve Model_hybrid:
    [Model_FSPsoln,Model_hybrid.fspOptions.bounds] = Model_hybrid.solve; 
    
    % Plot marginal distributions:
    Model_hybrid.makePlot(Model_FSPsoln,'marginals',[1:100:100],...
                       false,[1,2,3],{'linewidth',2})  
    Model_hybrid.makePlot(Model_FSPsoln,'margmovie',[],false,[101],...
                       {'linewidth',2},'movie.mp4',[],[],plotTitle='mRNA')