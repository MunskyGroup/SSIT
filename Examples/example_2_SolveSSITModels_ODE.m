%% SSIT/Examples/example_2_SolveSSITModels_ODE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Compute Ordinary Differential Equations (ODEs)
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
%% Ex(1): Use deterministic, ordinary differential equations (ODEs) 
% to average the time evolution of state space probabilities for 
% the bursting gene example model from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model:

    % Solve ODEs:
    Model = Model.solve(solver='ODE'); 

    % Plot ODE solutions:
    Model.plotODE(speciesNames=Model.species, timeVec=Model.tSpan,...
        lineProps={'linewidth',4}, AxisLabelSize=24, TickLabelSize=24,...
        Title='Busting Gene', TitleFontSize=32,...
        LegendFontSize=15, LegendLocation='southeast',...
        XLabel='Time', YLabel='Molecule Count')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use deterministic, ordinary differential equations (ODEs) 
% to average the time evolution of state space probabilities for 
% the time-varying STL1 yeast model from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1 Model:
    
    % Solve ODEs:
    STL1 = STL1.solve(solver='ODE'); 

    % Plot ODE solutions:
    STL1.plotODE(speciesNames=STL1.species, timeVec=STL1.tSpan,...
        lineProps={'linewidth',4}, Title='STL1', TitleFontSize=32,...
        AxisLabelSize=24, TickLabelSize=24,...
        LegendFontSize=15, LegendLocation='northwest',...
        XLabel='Time', YLabel='Molecule Count')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use deterministic, ordinary differential equations (ODEs) 
% to average the time evolution of state space probabilities for the
% 4-state time-varying STL1 yeast model from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4-state STL1 Model:

    % Set the ODE integrator (default 'ode23s'):
    STL1_4state.odeIntegrator = 'ode23s';
        
    % Solve ODEs:
    STL1_4state = STL1_4state.solve(solver='ODE'); 

    % Plot ODE solutions for mRNA:
    STL1_4state.plotODE(speciesNames=STL1_4state.species(5),...
        timeVec=STL1_4state.tSpan, lineProps={'linewidth',4},...
        TitleFontSize=26, Title='4-state STL1 (mRNA)',...
        AxisLabelSize=20, TickLabelSize=20, LegendFontSize=20,...
        LegendLocation='east', Colors=[0.23,0.67,0.20],...
        XLabel='Time', YLabel='Molecule Count')

    % Plot ODE solutions for the four gene states:
    STL1_4state.plotODE(speciesNames=STL1_4state.species(1:4),...
        timeVec=STL1_4state.tSpan, lineProps={'linewidth',4},...
        TitleFontSize=26, Title='4-state STL1 (gene states)',...
        AxisLabelSize=20, TickLabelSize=20, LegendFontSize=20,...
        LegendLocation='east', XLabel='Time', YLabel='Molecule Count')

    % Make a movie of the ODE solution being plotted:
    % STL1_4state.plotODE(makeMovie=true,...
    %     videoFileName='STL1_4state_ODE.mp4');


%% Save ODE models & solutions
saveNames = unique({'Model'
    'STL1'
    'STL1_4state'
    });
    
save('example_2_SolveSSITModels_ODE',saveNames{:})