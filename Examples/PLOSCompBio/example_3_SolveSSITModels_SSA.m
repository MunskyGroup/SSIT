%% example_3_SolveSSITModels_SSA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%% Compute Stochastic Simulation Algorithm (SSA) trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath(genpath('../../src'));

%% Preliminaries
% example_1_CreateSSITModels

% Load our models from example_1_CreateSSITModels and inspect them:
Model.summarizeModel
STL1.summarizeModel

% Set the times at which distributions will be computed:
Model.tSpan = linspace(0,20,200);
STL1.tSpan = linspace(0,20,200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Use Gillepsie's Stochastic Simulation Algorithm (SSA) 
% to solve the time evolution of state space probabilities for 
% the bursting gene example model from example_1_CreateSSITModels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    % 'nSimsPerExpt' is an SSA option that defaults to 100, sets the number
    % of simulations performed per experiment (set small number for demo)
    Model_SSA.ssaOptions.nSimsPerExpt=10;

    % 'verbose' defaults to false, prints completed sim number to screen
    Model_SSA.ssaOptions.verbose=true;
    
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

    %% Make a video of the SSA trajectories being plotted:
    makeSSAvideo(Model_SSAsoln, 'all', 100, Model_SSA.species, ...
        'Model_SSA_video')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use Gillepsie's Stochastic Simulation Algorithm (SSA) 
% to solve the time evolution of state space probabilities for 
% the time-varying STL1 yeast model from example_1_CreateSSITModels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1 Model:
    % Create a copy of the STL1 Model for SSAs:
    STL1_SSA = STL1;

    % Set solution scheme to SSA:
    STL1_SSA.solutionScheme = 'SSA';

    % 'nSimsPerExpt' is an SSA option that defaults to 100, sets the number
    % of simulations performed per experiment (set small number for demo)
    STL1_SSA.ssaOptions.nSimsPerExpt=10;

    % 'verbose' defaults to false, prints completed sim number to screen
    STL1_SSA.ssaOptions.verbose=true;
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_SSA'
    STL1_SSA = STL1_SSA.formPropensitiesGeneral('STL1_SSA');
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    STL1_SSA.tSpan = [-1,STL1_SSA.tSpan];

    % Set the initial time:
    STL1_SSA.initialTime = STL1_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1_SSAsoln = STL1_SSA.solve;
            
    % Plot SSA trajectories and means:
    plotSSA(STL1_SSAsoln, 'all', 100, STL1_SSA.species);

    %% Make a video of the SSA trajectories being plotted:
    makeSSAvideo(STL1_SSAsoln, 'all', 100, STL1_SSA.species, ...
        'STL1_SSA_video')