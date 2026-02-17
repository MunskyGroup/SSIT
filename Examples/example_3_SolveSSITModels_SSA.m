%% SSIT/Examples/example_3_SolveSSITModels_SSA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Compute Stochastic Simulation Algorithm (SSA) trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels
% clear
% close all
addpath(genpath('../src'));

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
%% Ex(1): Use Gillepsie's Stochastic Simulation Algorithm (SSA) 
% to solve the time evolution of state space probabilities for   
% the bursting gene example model from example_1_CreateSSITModels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run Gillepsie's Stochastic Simulation Algorithm (SSA) and analyse 
%% trajectories
%% Model:
    % Create a copy of the bursting gene model for SSA:
    Model_SSA = Model;

    % Set solution scheme to SSA:
    Model_SSA.solutionScheme = 'SSA';

    % 'nSimsPerExpt' is an SSA option that defaults to 100, sets the number
    % of simulations performed per experiment (set small number for demo)
    Model_SSA.ssaOptions.nSimsPerExpt=10;

    % 'verbose' defaults to false, prints completed sim number to screen
    Model_SSA.ssaOptions.verbose=true;
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    Model_SSA.tSpan = [-100,Model_SSA.tSpan];
    
    % Set the initial time:
    Model_SSA.initialTime = Model_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    Model_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    Model_SSA.Solutions = Model_SSA.solve;
    
    % Plot SSA trajectories and means:
    Model_SSA.plotSSA('all', 10, Model_SSA.species, {'linewidth',4}, ...
        Title="Bursting Gene", MeanOnly=true, TitleFontSize=32,...
        AxisLabelSize=24, TickLabelSize=24, LegendFontSize=20,...
        LegendLocation='east', XLabel='Time', YLabel='Molecule Count');

    %% Make a video of the SSA trajectories being plotted:
    % makeSSAvideo(Model_SSAsoln, 'all', 100, Model_SSA.species, ...
    %             'Model_SSA_video')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use Gillepsie's Stochastic Simulation Algorithm (SSA) 
% to solve the time evolution of state space probabilities for the 
% time-varying STL1 yeast model from example_1_CreateSSITModels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1 model:
    % Create a copy of the time-varying STL1 yeast model for SSA:
    STL1_SSA = STL1;

    % Set solution scheme to SSA:
    STL1_SSA.solutionScheme = 'SSA';

    % 'nSimsPerExpt' is an SSA option that defaults to 100, sets the number
    % of simulations performed per experiment (set small number for demo)
    STL1_SSA.ssaOptions.nSimsPerExpt=10;

    % 'verbose' defaults to false, prints completed sim number to screen
    STL1_SSA.ssaOptions.verbose=true;
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    STL1_SSA.tSpan = [-100,STL1_SSA.tSpan];

    % Set the initial time:
    STL1_SSA.initialTime = STL1_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1_SSA.Solutions = STL1_SSA.solve;
            
    % Plot SSA trajectories and means:
    STL1_SSA.plotSSA('all', 100, STL1_SSA.species, {'linewidth',4}, ...
        Title="STL1", MeanOnly=true, TitleFontSize=32,...
        AxisLabelSize=24, TickLabelSize=24,...
        LegendFontSize=20, LegendLocation='east',...
        XLabel='Time', YLabel='Molecule Count');

    %% Make a video of the SSA trajectories being plotted:
    % makeSSAvideo(STL1_SSAsoln, 'all', 100, STL1_SSA.species, ...
    %             'STL1_SSA_video')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use Gillepsie's Stochastic Simulation Algorithm (SSA) 
% to solve the time evolution of state space probabilities for the 4-state
% time-varying STL1 yeast model from example_1_CreateSSITModels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4-state STL1 model:
    % Create a copy of the time-varying STL1 yeast model for SSA:
    STL1_4state_SSA = STL1_4state;

    % Set solution scheme to SSA:
    STL1_4state_SSA.solutionScheme = 'SSA';

    % 'nSimsPerExpt' is an SSA option that defaults to 100, sets the number
    % of simulations performed per experiment (set small number for demo)
    STL1_4state_SSA.ssaOptions.nSimsPerExpt=10;

    % 'verbose' defaults to false, prints completed sim number to screen
    STL1_4state_SSA.ssaOptions.verbose=true;
       
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    STL1_4state_SSA.tSpan = [-100,STL1_4state_SSA.tSpan];

    % Set the initial time:
    STL1_4state_SSA.initialTime = STL1_4state_SSA.tSpan(1); 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1_4state_SSA.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1_4state_SSA.Solutions = STL1_4state_SSA.solve;
            
    % Plot SSA trajectories and means (mRNA):
    STL1_4state_SSA.plotSSA('all', 100, STL1_4state_SSA.species(5),...
        {'linewidth',4}, Title="4-state STL1 (mRNA)", MeanOnly=false,...
        TitleFontSize=26, AxisLabelSize=20, TickLabelSize=20,...
        LegendFontSize=20, LegendLocation='northeast', HistTime=20,...
        XLabel='Time', YLabel='Molecule Count', Colors=[0.23,0.67,0.20]);  

    % Plot SSA trajectories and means (gene states):
    STL1_4state_SSA.plotSSA('all', 100, STL1_4state_SSA.species(1:4),...
        {'linewidth',4}, Title="4-state STL1 (gene states)",...
        MeanOnly=true, TitleFontSize=26, AxisLabelSize=20,...
        TickLabelSize=20, LegendFontSize=20, LegendLocation='east',...
        XLabel='Time', YLabel='Molecule Count'); 

    %% Make a video of the SSA trajectories being plotted:
    % makeSSAvideo(STL1_4state_SSAsoln, 'all', 100, ...
    %              STL1_4state_SSA.species, 'STL1_SSA_video_4state')

%% Save SSA models & solutions
saveNames = unique({'Model_SSA'
    'STL1_SSA'
    'STL1_4state_SSA'
    });
    
save('example_3_SolveSSITModels_SSA',saveNames{:})