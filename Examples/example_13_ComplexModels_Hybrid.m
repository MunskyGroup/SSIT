%% example_13_SolveSSITModels_Hybrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Solve a model by using a hybrid deterministic-stochastic approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels
%clear
%close all
addpath(genpath('../src'));

% example_1_CreateSSITModels 

% View model summaries:
STL1.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust a model to treat some species (i.e., upstream reactions) use an 
% ODE formulation, while having other species (i.e., downstream species) 
% evolving in a discrete stochastic manner. This runs significantly faster 
% than full FSP solutions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create copies of our models:
STL1_hybrid = STL1;

% Set the times at distributions will be computed:
STL1_hybrid.tSpan = linspace(0,20,200);

% Set 'useHybrid' to true:
STL1_hybrid.useHybrid = true;

% Define which species will be solved by ODEs:
STL1_hybrid.hybridOptions.upstreamODEs = {'offGene','onGene'};

% Set solution scheme to FSP:
STL1_hybrid.solutionScheme = 'FSP';

% View summary of hybrid STL1 model, expecting:
    %   Species:
    %       offGene; IC = 1;  upstream ODE
    %       onGene; IC = 0;  upstream ODE
    %       mRNA; IC = 0;  discrete stochastic
STL1_hybrid.summarizeModel

% Optionally, define a custom constraint on the species: 
% (e.g.,'offGene'+'onGene')
STL1_hybrid.customConstraintFuns = [];

% Set FSP 1-norm error tolerance:
STL1_hybrid.fspOptions.fspTol = 1e-4; 
    
% Guess initial bounds on FSP StateSpace:
STL1_hybrid.fspOptions.bounds = [2,2,400];

% This function compiles and stores the given reaction propensities  
% into symbolic expression functions that use sparse matrices to  
% operate on the system based on the current state. The functions are 
% stored with the given prefix, in this case, 'STL1_hybrid' and
% 'STL1_4state_hybrid'
STL1_hybrid = STL1_hybrid.formPropensitiesGeneral('STL1_hybrid');

% Have FSP approximate the steady state for the initial distribution 
% by finding the eigenvector corresponding to the smallest magnitude 
% eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
STL1_hybrid.fspOptions.initApproxSS = false; 
    
% Solve STL1_hybrid:
[STL1_hybrid_FSPsoln,STL1_hybrid.fspOptions.bounds] = STL1_hybrid.solve; 
    
% Plot marginal distributions:
STL1_hybrid.makePlot(STL1_hybrid_FSPsoln,'meansAndDevs')  
STL1_hybrid.makePlot(STL1_hybrid_FSPsoln,'margmovie',[],false,[101],...
                  {'linewidth',2},'STL1_hybrid.mp4',[],[],plotTitle='mRNA')