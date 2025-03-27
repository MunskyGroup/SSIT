%% example_2_SolveSSITModels
% Example script to show how to solve the time evolution of state space 
% probabilities for a reaction system where processes are considered:  
% * Deterministic, using ordinary differential equations (ODEs) to average;
% * Stochastic, using stochastic simulation algorithm (SSA) trajectories 
% * or Finite State Projection (FSP) of the Chemical Master Equation (CME);
% * Hybrid, using upstream ODEs and downstream FSP.
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

%% Ex.(1) Compute Ordinary Differential Equations (ODEs)
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
    ODEsoln = Model_ODE.solve; 
    plotODE(ODEsoln,Model_ODE.species)

%% STL1 Model:
    % Create a copy of the STL1 Model for ODEs:
    STL1Model_ODE = STL1Model;

    % Set solution scheme to 'ode':
    STL1Model_ODE.solutionScheme = 'ode';
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1Model_ODE'
    STL1Model_ODE = STL1Model_ODE.formPropensitiesGeneral('STL1Model_ODE');
    
    % Solve ODE and make plots:
    STL1_ODEsoln = STL1Model_ODE.solve; 
    plotODE(STL1_ODEsoln,STL1Model_ODE.species)

%% Ex.(2): Run Gillepsie's Stochastic Simulation Algorithm (SSA) and
%%         analyse trajectories
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
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). This can cause long run times.
    Model_SSA.tSpan = [-10,Model_SSA.tSpan];
    
    % Set the initial time:
    Model_SSA.initialTime = Model_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    Model_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    SSAsoln = Model_SSA.solve;
            
    % Plot SSA results:
    plotSSA(SSAsoln, 'all', 1200, Model_SSA.species);

%% STL1 Model:
    % Create a copy of the STL1 Model for SSAs:
    STL1Model_SSA = STL1Model;

    % Set solution scheme to SSA:
    STL1Model_SSA.solutionScheme = 'SSA';
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1Model_SSA'
    STL1Model_SSA = STL1Model_SSA.formPropensitiesGeneral('STL1Model_SSA');
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). This can cause long run times.
    STL1Model_SSA.tSpan = [-10,STL1Model_SSA.tSpan];

    % Set the initial time:
    STL1Model_SSA.initialTime = STL1Model_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1Model_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1_SSAsoln = STL1Model_SSA.solve;
            
    % Plot SSA tesults:
    plotSSA(STL1_SSAsoln, 'all', 1200, STL1Model_SSA.species);


%% Ex.(3) Use the Finite State Projection (FSP) approximation of the CME
%% Model:
    % Create a copy of the Model for FSP:
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
    [FSPsoln,Model_FSP.fspOptions.bounds] = Model_FSP.solve; 
    
    % Plot marginal distributions:
    Model_FSP.makePlot(FSPsoln,'marginals',[1:100:100],...
                       false,[1,2,3],{'linewidth',2})  
    Model_FSP.makePlot(FSPsoln,'margmovie',[],false,[101],...
                       {'linewidth',2},'movie.mp4',[1,1,0.6],[2,3])  

%% STL1Model:
    % Create a copy of the STL1 Model for FSP:
    STL1Model_FSP = STL1Model;
    
    % Ensure the solution scheme is set to FSP (default):
    STL1Model_FSP.solutionScheme = 'FSP';  

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1Model_FSP'
    STL1Model_FSP = STL1Model_FSP.formPropensitiesGeneral('STL1Model_FSP');
    
    % Set FSP 1-norm error tolerance:
    STL1Model_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1Model_FSP.fspOptions.bounds = [2,2,400];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1Model_FSP.fspOptions.initApproxSS = true; 
    
    % Solve Model:
    [STL1_FSPsoln,STL1Model_FSP.fspOptions.bounds] = STL1Model_FSP.solve; 
    
    % Plot marginal distributions:
    STL1Model_FSP.makePlot(STL1_FSPsoln,'marginals',[1:100:100],...
                           false,[1,2,3],{'linewidth',2})  
    STL1Model_FSP.makePlot(STL1_FSPsoln,'margmovie',[],false,[101],...
                           {'linewidth',2},'movie.mp4',[1,1,1],[2,3])  