%% SSIT/Examples/example_3_SolveSSITModels_SSA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.2.3: Finding and visualizing master equation solutions: SSA
%   * Compute Stochastic Simulation Algorithm (SSA) trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels
% clear
% close all

% example_1_CreateSSITModels

% Load the models created in example_1_CreateSSITModels
load('example_1_CreateSSITModels.mat')

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
    % Set the number of simulations:
    Model.ssaOptions.nSims=1000;

    % 'verbose' defaults to false, prints completed sim number to screen.
    Model.ssaOptions.verbose=true;
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    Model.tSpan = [-100,Model.tSpan];
    
    % Set the initial time:
    Model.initialTime = -100; 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    Model.ssaOptions.useParallel = true;
    
    % Run SSA:
    Model = Model.solve(solver='SSA');
    
    % Plot SSA trajectories and means:
    Model.plotSSA(speciesIdx='all', numTraj=10,...
        speciesNames=Model.species, lineProps={'linewidth',4}, ...
        Title="Bursting Gene", MeanOnly=true, TitleFontSize=32,...
        AxisLabelSize=24, TickLabelSize=24, LegendFontSize=20,...
        LegendLocation='east', XLabel='Time', YLabel='Molecule Count');

    %% Make a video of the SSA trajectories being plotted:
    % Model.plotSSA(makeMovie=true,
    %             'Model_video')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use Gillepsie's Stochastic Simulation Algorithm (SSA) 
% to solve the time evolution of state space probabilities for the 
% time-varying STL1 yeast model from example_1_CreateSSITModels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1 model:
    % Set the number of simulations:
    STL1.ssaOptions.nSims=1000;
    
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    STL1.tSpan = [-100,STL1.tSpan];

    % Set the initial time:
    STL1.initialTime = -100; 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1 = STL1.solve(solver='SSA');
            
    % Plot SSA trajectories and means:
    STL1.plotSSA(speciesIdx='all', numTraj=100,...
        speciesNames=STL1.species, lineProps={'linewidth',4},...
        Title="STL1", MeanOnly=true, TitleFontSize=32,...
        AxisLabelSize=24, TickLabelSize=24,...
        LegendFontSize=20, LegendLocation='east',...
        XLabel='Time', YLabel='Molecule Count');

    %% Make a video of the SSA trajectories being plotted:
    % makeSSAvideo(STL1.Solutions, 'all', 100, STL1.species, ...
    %             'STL1_video')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use Gillepsie's Stochastic Simulation Algorithm (SSA) 
% to solve the time evolution of state space probabilities for the 4-state
% time-varying STL1 yeast model from example_1_CreateSSITModels  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4-state STL1 model:
    % Set the number of simulations:
    STL1_4state.ssaOptions.nSims=1000;
       
    % A negative initial time is used to allow model to equilibrate 
    % before starting (burn-in). Large burn-in times cause long run times.
    STL1_4state.tSpan = [-100,STL1_4state.tSpan];

    % Set the initial time:
    STL1_4state.initialTime = -100; 
    
    % Run iterations in parallel with multiple cores, or execute serially:
    STL1_4state.ssaOptions.useParallel = true;
    
    % Run SSA:
    STL1_4state = STL1_4state.solve(solver='SSA');
            
    % Plot SSA trajectories and means (mRNA):
    STL1_4state.plotSSA(speciesIdx='all', numTraj=100,...
        speciesNames=STL1_4state.species(5), MeanOnly=false,...
        lineProps={'linewidth',4}, Title="4-state STL1 (mRNA)",...
        TitleFontSize=26, AxisLabelSize=20, TickLabelSize=20,...
        LegendFontSize=20, LegendLocation='northeast', HistTime=20,...
        XLabel='Time', YLabel='Molecule Count', Colors=[0.23,0.67,0.20]);  

    % Plot SSA trajectories and means (gene states):
    STL1_4state.plotSSA(speciesIdx='all', numTraj=100,...
        speciesNames=STL1_4state.species(1:4),...
        lineProps={'linewidth',4}, Title="4-state STL1 (gene states)",...
        MeanOnly=true, TitleFontSize=26, AxisLabelSize=20,...
        TickLabelSize=20, LegendFontSize=20, LegendLocation='east',...
        XLabel='Time', YLabel='Molecule Count', makeMovie=false); 

    %% Make a video of the SSA trajectories being plotted:
    % STL1_4state.plotSSA(makeMovie=true,videoFileName='STL1_video_4state',...
    %  speciesIdx='all',numTraj=100,speciesNames=STL1_4state.species(1:5))

%% Save SSA models & solutions
saveNames = unique({'Model'
    'STL1'
    'STL1_4state'
    });
    
save('example_3_SolveSSITModels_SSA',saveNames{:})