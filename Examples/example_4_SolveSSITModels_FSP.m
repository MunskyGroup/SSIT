%% example_4_SolveSSITModels_FSP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Compute Finite State Projection (FSP) solutions
%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels
% clear
% close all

% example_1_CreateSSITModels

% Load the models created in example_1_CreateSSITModels
% load('example_1_CreateSSITModels.mat')

% View model summaries:
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

% Set the times at which distributions will be computed:
Model.tSpan = linspace(0,50,101);
STL1.tSpan = linspace(0,50,101);
STL1_4state.tSpan = linspace(0,50,101);

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
    
    % Set FSP 1-norm error tolerance:
    Model_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    Model_FSP.fspOptions.bounds = [0,0,0,1,1,200];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    Model_FSP.fspOptions.initApproxSS = false; 

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'Model_FSP':
    Model_FSP = Model_FSP.formPropensitiesGeneral('Model_FSP');
    
    % Solve with FSP:
    [~,~,Model_FSP] = Model_FSP.solve; 
    
    % Plot marginal distributions at t=20:  
    Model_FSP.plotFSP(plotType='meansAndDevs',Title='Bursting Gene (FSP)')
                       
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

    % Set FSP 1-norm error tolerance:
    STL1_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1_FSP.fspOptions.bounds = [0,0,0,1,1,200];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1_FSP.fspOptions.initApproxSS = false; 

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_FSP':
    STL1_FSP = STL1_FSP.formPropensitiesGeneral('STL1_FSP');
    
    % Solve with FSP:
    [~,~,STL1_FSP] = STL1_FSP.solve; 
    
    % Plot marginal distributions at t=20:  
    STL1_FSP.plotFSP(plotType='meansAndDevs',Title='STL1 (FSP)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use the stochastic Finite State Projection (FSP) 
% approximation of the Chemical Master Equation (CME) to solve the time 
% evolution of state space probabilities for the 4-state 
% time-varying STL1 yeast model from example_1_CreateSSITModels 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1 (4-state):
    % Create a copy of the time-varying STL1 yeast model for FSP:
    STL1_4state_FSP = STL1_4state;
    
    % Ensure the solution scheme is set to FSP (default):
    STL1_4state_FSP.solutionScheme = 'FSP';  
    
    % Set FSP 1-norm error tolerance:
    STL1_4state_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1_4state_FSP.fspOptions.bounds = [1,1,1,1,200];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1_4state_FSP.fspOptions.initApproxSS = true; 

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_4state_FSP':
    STL1_4state_FSP = ...
        STL1_4state_FSP.formPropensitiesGeneral('STL1_4state_FSP');
    
    % Solve Model:
    [~,~,STL1_4state_FSP] = STL1_4state_FSP.solve; 
    
    %% Plots for FSP solutions:
    % Means only:
    % STL1_4state_FSP.plotFSP(STL1_4state_FSPsoln,...
    %     STL1_4state_FSP.species, 'means')

    % Means and standard deviations:
    STL1_4state_FSP.plotFSP(speciesNames=STL1_4state_FSP.species(5),...
        plotType='meansAndDevs', lineProps={'linewidth',4},...
        Title='4-state STL1 (mRNA)', TitleFontSize=26,...
        Colors=[0.23,0.67,0.2], AxisLabelSize=20, TickLabelSize=20,...
        XLabel='Time', YLabel='Molecule Count',...
        LegendFontSize=20, LegendLocation='northeast');

    % Marginal distributions:
    STL1_4state_FSP.plotFSP(speciesNames=STL1_4state_FSP.species(5),...
        plotType='marginals', indTimes=[1,12,24,50,101],...
        lineProps={'linewidth',3}, Colors=[0.23,0.67,0.2], XLim=[0,100])

    % Joint distributions (warning: can be slow for many parameters!):
    % STL1_4state_FSP.plotFSP(STL1_4state_FSPsoln,...
    %     STL1_4state_FSP.species(5), 'joints')


%% Save FSP models & solutions
saveNames = unique({
    'Model_FSP'
    'STL1_FSP'
    'STL1_4state_FSP'
    });
    
save('example_4_SolveSSITModels_FSP',saveNames{:})