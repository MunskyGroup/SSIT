%% example_4_SolveSSITModels_FSP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Compute Finite State Projection (FSP) solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels
%clear
%close all
%addpath(genpath('../../'));

% example_1_CreateSSITModels

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
%% Ex(1): Use the stochastic Finite State Projection (FSP) 
% approximation of the Chemical Master Equation (CME) to solve the time 
% evolution of state space probabilities for the bursting gene example 
% model from example_1_CreateSSITModels 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model:
    % Create a copy of the bursting gene model for FSP:
    Model_FSP = Model;
    
    % Ensure the solution scheme is set to FSP (default):
    Model_FSP.solutionScheme = 'FSP';  

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'Model_FSP'
    Model_FSP = Model_FSP.formPropensitiesGeneral('Model_FSP');
    
    % Set FSP 1-norm error tolerance:
    Model_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    Model_FSP.fspOptions.bounds = [2,2,400];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    Model_FSP.fspOptions.initApproxSS = false; 
    
    % Solve with FSP:
    [Model_FSPsoln,Model_FSP.fspOptions.bounds] = Model_FSP.solve; 
    
    % Plot marginal distributions:
    Model_FSP.makePlot(Model_FSPsoln,'marginals',[1:100:100],...
                       false,[1,2,3],{'linewidth',2})  
    Model_FSP.makePlot(Model_FSPsoln,'meansAndDevs') 
    Model_FSP.makePlot(Model_FSPsoln,'margmovie',[],false,[101],...
                        {'linewidth',2},'Model_FSP.mp4',[1,1,0.5],[2,3]) 
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use the stochastic Finite State Projection (FSP) 
% approximation of the Chemical Master Equation (CME) to solve the time 
% evolution of state space probabilities for the time-varying STL1 yeast 
% model from example_1_CreateSSITModels 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1:
    % Create a copy of the time-varying STL1 yeast model for FSP:
    STL1_FSP = STL1;
    
    % Ensure the solution scheme is set to FSP (default):
    STL1_FSP.solutionScheme = 'FSP';  

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_FSP'
    STL1_FSP = STL1_FSP.formPropensitiesGeneral('STL1_FSP');
    
    % Set FSP 1-norm error tolerance:
    STL1_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1_FSP.fspOptions.bounds = [2,2,400];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1_FSP.fspOptions.initApproxSS = false; 
    
    % Solve Model:
    [STL1_FSPsoln,STL1_FSP.fspOptions.bounds] = STL1_FSP.solve; 
    
    % Plot marginal distributions:
    STL1_FSP.makePlot(STL1_FSPsoln,'marginals',[1:100:100],...
                           false,[1,2,3],{'linewidth',2}) 
    STL1_FSP.makePlot(STL1_FSPsoln,'meansAndDevs')  
    STL1_FSP.makePlot(STL1_FSPsoln,'margmovie',[],false,[101],...
                           {'linewidth',2},'STL1_FSP.mp4',[1,1,0.5],[2,3])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use the stochastic Finite State Projection (FSP) 
% approximation of the Chemical Master Equation (CME) to solve the time 
% evolution of state space probabilities for the 4-state 
% time-varying STL1 yeast model from example_1_CreateSSITModels 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1:
    % Create a copy of the time-varying STL1 yeast model for FSP:
    STL1_4state_FSP = STL1_4state;
    
    % Ensure the solution scheme is set to FSP (default):
    STL1_4state_FSP.solutionScheme = 'FSP';  

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_FSP'
    STL1_4state_FSP = ...
        STL1_4state_FSP.formPropensitiesGeneral('STL1_4state_FSP');
    
    % Set FSP 1-norm error tolerance:
    STL1_4state_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1_4state_FSP.fspOptions.bounds = [2,2,2,2,400];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1_4state_FSP.fspOptions.initApproxSS = false; 
    
    % Solve Model:
    [STL1_4state_FSPsoln,STL1_4state_FSP.fspOptions.bounds] = ...
        STL1_4state_FSP.solve; 
    
    % Plot marginal distributions:
    STL1_4state_FSP.makePlot(STL1_4state_FSPsoln,'marginals',...
                             [1:100:100],false,[1,2,3,4,5],{'linewidth',2})
    STL1_4state_FSP.makePlot(STL1_4state_FSPsoln,'meansAndDevs')  
    STL1_4state_FSP.makePlot(STL1_4state_FSPsoln,'margmovie',[],false,...
                             [101],{'linewidth',2},'STL1_4state_FSP.mp4')

%% Save FSP models & solutions
saveNames = unique({'Model_FSP'
    'Model_FSPsoln'
    'STL1_FSP'
    'STL1_FSPsoln'
    'STL1_4state_FSP'
    'STL1_4state_FSPsoln'
    });
    
save('example_4_SolveSSITModels_FSP',saveNames{:})