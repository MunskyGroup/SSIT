%% SSIT/Examples/example_13_SolveSSITModels_Hybrid

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
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust a model to treat some species (i.e., upstream reactions) using an 
% ODE formulation, while having other species (i.e., downstream species) 
% evolve in a discrete stochastic manner using FSP. This runs significantly 
% faster than full FSP solutions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a copy of our models:
STL1_hybrid = STL1_4state;

% Set the times at distributions will be computed:
STL1_hybrid.tSpan = linspace(0,50,200);

% Set 'useHybrid' to true:
STL1_hybrid.useHybrid = true;

% Define which species will be solved by ODEs:
STL1_hybrid.hybridOptions.upstreamODEs = {'g1','g2','g3','g4'};

% Set solution scheme to FSP:
STL1_hybrid.solutionScheme = 'FSP';

% View summary of hybrid STL1 model, expecting:
    %   Species:
    %       g1; IC = 1;  upstream ODE
    %       g2; IC = 0;  upstream ODE
    %       g3; IC = 0;  upstream ODE
    %       g4; IC = 0;  upstream ODE
    %       mRNA; IC = 0;  discrete stochastic
STL1_hybrid.summarizeModel

% Optionally, define a custom constraint on the species: 
% (e.g.,'offGene'+'onGene')
STL1_hybrid.customConstraintFuns = [];

% Set FSP 1-norm error tolerance:
STL1_hybrid.fspOptions.fspTol = 1e-4; 
    
% Guess initial bounds on FSP StateSpace:
STL1_hybrid.fspOptions.bounds = [1,1,1,1,300];

% This function compiles and stores the given reaction propensities  
% into symbolic expression functions that use sparse matrices:
STL1_hybrid = STL1_hybrid.formPropensitiesGeneral('STL1_hybrid');

% Have FSP approximate the steady state for the initial distribution:
STL1_hybrid.fspOptions.initApproxSS = true; 
    
% Solve STL1_hybrid:
[STL1_hybrid_FSPsoln,STL1_hybrid.fspOptions.bounds] = STL1_hybrid.solve; 

%% Plots for FSP solutions:
    % Means and standard deviations:
    STL1_hybrid.plotFSP(STL1_hybrid_FSPsoln,...
        STL1_hybrid.species(5), 'meansAndDevs', [], [],...
        {'linewidth',4}, Title='4-state STL1 (mRNA)', TitleFontSize=24,...
        Colors=[0.23,0.67,0.2], AxisLabelSize=18, TickLabelSize=18,...
        XLabel='Time', YLabel='Molecule Count',...
        LegendFontSize=15, LegendLocation='northeast');

    % Marginal distributions:
    STL1_hybrid.plotFSP(STL1_hybrid_FSPsoln,...
        STL1_hybrid.species(5), 'marginals', [1,12,24,50,101,200],...
        [], {'linewidth',3}, Colors=[0.23,0.67,0.2], AxisLabelSize=18,...
        TickLabelSize=18, XLim=[0,100])