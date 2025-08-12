%% example_2_SolveSSITModels_ODE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Compute Ordinary Differential Equations (ODEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels
%clear
%close all
addpath(genpath('../../'));

%example_1_CreateSSITModels

% Load the models created in example_1_CreateSSITModels
load('example_1_CreateSSITModels.mat')

% View model summaries:
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

% Set the times at which distributions will be computed:
Model.tSpan = linspace(0,20,200);
STL1.tSpan = linspace(0,20,200);
STL1_4state.tSpan = linspace(0,20,200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Use deterministic, ordinary differential equations (ODEs) 
% to average the time evolution of state space probabilities for 
% the bursting gene example model from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model:
    % Create a copy of the bursting gene model for ODEs:
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

    %% Make a movie of the ODE solution being plotted:
    makeODEmovie(Model_ODEsoln, Model_ODE.species, Model_ODE.tSpan, ...
                'Model_ODE.mp4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use deterministic, ordinary differential equations (ODEs) 
% to average the time evolution of state space probabilities for 
% the time-varying STL1 yeast model from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1 Model:
    % Create a copy of the time-varying STL1 yeast model for ODEs:
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

    %% Make a movie of the ODE solution being plotted:
    makeODEmovie(STL1_ODEsoln, STL1_ODE.species, STL1_ODE.tSpan, ...
                'STL1_ODE.mp4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use deterministic, ordinary differential equations (ODEs) 
% to average the time evolution of state space probabilities for the
% 4-state time-varying STL1 yeast model from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1 Model:
    % Create a copy of the time-varying STL1 yeast model for ODEs:
    STL1_4state_ODE = STL1_4state;

    % Set solution scheme to 'ode':
    STL1_4state_ODE.solutionScheme = 'ode';
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_ODE'
    STL1_4state_ODE = ...
        STL1_4state_ODE.formPropensitiesGeneral('STL1_4state_ODE');
    
    % Solve ODE and make plots:
    STL1_4state_ODEsoln = STL1_4state_ODE.solve; 
    plotODE(STL1_4state_ODEsoln,STL1_4state_ODE.species,...
            STL1_4state_ODE.tSpan)

    %% Make a movie of the ODE solution being plotted:
    makeODEmovie(STL1_4state_ODEsoln, STL1_4state_ODE.species, ...
                 STL1_4state_ODE.tSpan, 'STL1_4state_ODE.mp4');


%% Save ODE models & solutions
saveNames = unique({'Model_ODE'
    'STL1_ODE'
    'STL1_4state_ODE'
    'Model_ODEsoln'
    'STL1_ODEsoln'
    'STL1_4state_ODEsoln'
    });
    
save('example_2_SolveSSITModels_ODE',saveNames{:})