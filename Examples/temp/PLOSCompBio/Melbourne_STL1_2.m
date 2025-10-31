%% example_1_CreateSSITModels
%  Example script to show how to create models in SSIT.
clear
close all
addpath(genpath('../../src'));

example_1a_CreateSSITModels

%% Ex.(2) Create an SSIT model for real, time-varying experimental data: 
%%           yeast STL1_2 data provided by the Vanderbilt Q_BIO Group

% Create, solve, and fit a CME model to single-cell smFISH data. 
% The dataset is from yeast STL1_2 collected in Dr. Gregor Neuert's 
% laboratory at Vanderbilt University.

% Copy the Bursting Gene Model above to a new model, 'STL1_2'.  
% The species, reactions, and stochiometries remain the same; 
% however, a change is made to the propensity for the first reaction 
% (gene activation) so that its rate is controlled by a time-varying 
% MAPK input signal. As a result, new parameters must also be added.

% Create a copy of the simple Bursting Gene model from above:
STL1_2 = STL1;


% Define the time-varying TF/MAPK input signal:
STL1_2.inputExpressions = {'IHog',...
                         '(a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0))';...
                         'IHog',...
                         '(1000+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>5))'};

% Print a summary of STL1_2 Model:
STL1_2.summarizeModel
STL1_2.tSpan = linspace(0,20,200);

%% STL1_2 Model:
    % Create a copy of the STL1_2 Model for ODEs:
    STL1_2_ODE = STL1_2;

    % Set solution scheme to 'ode':
    STL1_2_ODE.solutionScheme = 'ode';
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_2_ODE'
    STL1_2_ODE = STL1_2_ODE.formPropensitiesGeneral('STL1_2_ODE');
    
    % Solve ODE and make plots:
    STL1_2_ODEsoln = STL1_2_ODE.solve; 
    plotODE(STL1_2_ODEsoln,STL1_2_ODE.species,STL1_2_ODE.tSpan)

    %% Make a movie of the ODE solution being plotted:
    makeODEmovie(STL1_2_ODEsoln, STL1_2_ODE.species, STL1_2_ODE.tSpan, ...
                'STL1_2_ODE.mp4');

%% Set solution scheme to SSA:
    STL1_SSA_2 = STL1_2; 

    STL1_SSA_2.solutionScheme = 'SSA';

    % 'nSimsPerExpt' is an SSA option that defaults to 100, sets the number
    % of simulations performed per experiment (set small number for demo)
    STL1_SSA_2.ssaOptions.nSimsPerExpt=10;

    % 'verbose' defaults to false, prints completed sim number to screen
    STL1_SSA_2.ssaOptions.verbose=true;
    
    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_SSA_2'
    STL1_SSA_2 = STL1_SSA_2.formPropensitiesGeneral('STL1_SSA_2');
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    STL1_SSA_2.tSpan = [-1,STL1_SSA_2.tSpan];

    % Set the initial time:
    STL1_SSA_2.initialTime = STL1_SSA_2.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1_SSA_2.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1_SSA_2soln = STL1_SSA_2.solve;
            
    % Plot SSA trajectories and means:
    plotSSA(STL1_SSA_2soln, 'all', 100, STL1_SSA_2.species);

    %% Make a video of the SSA trajectories being plotted:
    makeSSAvideo(STL1_SSA_2soln, 'all', 100, STL1_SSA_2.species, ...
        'STL1_SSA_2_video')

    %% STL1:
    % Create a copy of the STL1 Model for FSP:
    STL1_FSP_2 = STL1_2;
    
    % Ensure the solution scheme is set to FSP (default):
    STL1_FSP_2.solutionScheme = 'FSP';  

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_FSP_2'
    STL1_FSP_2 = STL1_FSP_2.formPropensitiesGeneral('STL1_FSP_2');
    
    % Set FSP 1-norm error tolerance:
    STL1_FSP_2.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1_FSP_2.fspOptions.bounds = [2,2,400];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1_FSP_2.fspOptions.initApproxSS = false; 
    
    % Solve Model:
    [STL1_FSP_2soln,STL1_FSP_2.fspOptions.bounds] = STL1_FSP_2.solve; 
    
    % Plot marginal distributions:
    STL1_FSP_2.makePlot(STL1_FSP_2soln,'marginals',[1:100:100],...
                           false,[1,2,3],{'linewidth',2})  
    STL1_FSP_2.makePlot(STL1_FSP_2soln,'margmovie',[],false,[101],...
                           {'linewidth',2},'STL1_FSP_2.mp4',[1,1,0.5],[2,3]) 
