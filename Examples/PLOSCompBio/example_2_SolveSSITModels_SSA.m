%% example_2_SolveSSITModels_SSA
% Example script to show how to solve the time evolution of state space 
% probabilities for a reaction system where processes are considered:  
% * Stochastic, using stochastic simulation algorithm (SSA) trajectories.
clear
close all
addpath(genpath('../../'));

%% Preliminaries
% Load our models from example_1_CreateSSITModels and inspect them:
example_1_CreateSSITModels
Model.summarizeModel
STL1Model.summarizeModel

% Set the times at distributions will be computed:
Model.tSpan = linspace(0,20,200);
STL1Model.tSpan = linspace(0,20,200);

%% Run Gillepsie's Stochastic Simulation Algorithm (SSA) and analyse 
%% trajectories
%% Model:
    % Create a copy of the Model for SSAs:
    Model_SSA = Model;

    % Set solution scheme to SSA:
    Model_SSA.solutionScheme = 'SSA';
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'Model_SSA'
    Model_SSA = Model_SSA.formPropensitiesGeneral('Model_SSA');

    Model_SSA.ssaOptions.verbose=true;
    Model_SSA.ssaOptions.nSimsPerExpt=10;

    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    Model_SSA.tSpan = [-1,Model_SSA.tSpan];
    
    % Set the initial time:
    Model_SSA.initialTime = Model_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    Model_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    Model_SSAsoln = Model_SSA.solve;
            
    % Plot SSA trajectories and means:
    plotSSA(Model_SSAsoln, 'all', 100, Model_SSA.species);

%% STL1 Model:
    % Create a copy of the STL1 Model for SSAs:
    STL1Model_SSA = STL1Model;

    % Set solution scheme to SSA:
    STL1Model_SSA.solutionScheme = 'SSA';

    STL1Model_SSA.ssaOptions.verbose=true;
    STL1Model_SSA.ssaOptions.nSimsPerExpt=10;
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1Model_SSA'
    STL1Model_SSA = STL1Model_SSA.formPropensitiesGeneral('STL1Model_SSA');
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    STL1Model_SSA.tSpan = [-1,STL1Model_SSA.tSpan];

    % Set the initial time:
    STL1Model_SSA.initialTime = STL1Model_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1Model_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1_SSAsoln = STL1Model_SSA.solve;
            
    % Plot SSA trajectories and means:
    plotSSA(STL1_SSAsoln, 'all', 2100, STL1Model_SSA.species);