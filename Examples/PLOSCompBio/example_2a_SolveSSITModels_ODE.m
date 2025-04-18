%% example_2_SolveSSITModels_ODE
% Example script to show how to solve the time evolution of state space 
% probabilities for a reaction system where processes are considered:  
% * Deterministic, using ordinary differential equations (ODEs) to average.
clear
close all
addpath(genpath('../../src'));

%% Preliminaries
% Load our models from example_1_CreateSSITModels and inspect them:
example_1a_CreateSSITModels
Model.summarizeModel
STL1.summarizeModel

% Set the times at distributions will be computed:
Model.tSpan = linspace(0,20,200);
STL1.tSpan = linspace(0,20,200);

%% Compute Ordinary Differential Equations (ODEs)
%% Model:
    % Create a copy of the Model for ODEs:
    Model_ODE = Model;

    % Set solution scheme to 'ode':
    Model_ODE.solutionScheme = 'ode';
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'Model_ODE'
    Model_ODE = Model_ODE.formPropensitiesGeneral('Model_ODE');
    
    % Solve ODE and make plots:
    Model_ODEsoln = Model_ODE.solve; 
    plotODE(Model_ODEsoln,Model_ODE.species,Model_ODE.tSpan)

%% STL1 Model:
    % Create a copy of the STL1 Model for ODEs:
    STL1_ODE = STL1;

    % Set solution scheme to 'ode':
    STL1_ODE.solutionScheme = 'ode';
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_ODE'
    STL1_ODE = STL1_ODE.formPropensitiesGeneral('STL1_ODE');
    
    % Solve ODE and make plots:
    STL1_ODEsoln = STL1_ODE.solve; 
    plotODE(STL1_ODEsoln,STL1_ODE.species,STL1_ODE.tSpan)