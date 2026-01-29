%% SSIT/Examples/example_5_SolveSSITModels_EscapeTimes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Solve a first-passage time problem (escape times)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels
% clear
% close all
addpath(genpath('../'));

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
%% Ex(1): Solve escape times for the bursting gene example model 
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model:
    % Create a copy of the bursting gene model:
    Model_escape = Model;
    
    %% Specify a boundary for the escape calculation
    % Calculate the time until the mRNA concentration reaches 5
    Model_escape.fspOptions.escapeSinks.f = {'mRNA'};
    Model_escape.fspOptions.verbose = false;
    Model_escape.fspOptions.escapeSinks.b = 5;
    Model_escape = Model_escape.formPropensitiesGeneral('Model_escape');
    [~,~,Model_escape] = Model_escape.solve;

    % Plot the CDF and PDF
    Model_escape.plotFSP(Model_escape.Solutions, [],...
        "escapeTimes", [], [], {'linewidth',3},...
        Title="Bursting Gene (mRNA)", Colors=[0.93,0.69,0.13],...
        LegendLocation="southeast");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve escape times for the time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1:
    % Create a copy of the time-varying STL1 yeast model:
    STL1_escape = STL1;
    
    % Solve for the escape time:
    STL1_escape.fspOptions.escapeSinks.f = {'offGene','onGene'};
    STL1_escape.fspOptions.escapeSinks.b = [0.5;0.5];
    STL1_escape = STL1_escape.formPropensitiesGeneral('STL1_escape');
    [~,~,STL1_escape] = STL1_escape.solve;

    % Plot the CDF and PDF
    STL1_escape.plotFSP(STL1_escape.Solutions,[],"escapeTimes",[],[],...
        {'linewidth',3}, Title="STL1 (offGene, onGene)",...
        LegendLocation="east");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve escape times for the 4-state time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4-state STL1:
    % Create a copy of the time-varying STL1 yeast model:
    STL1_4state_escape = STL1_4state;


    % This should not be required, since propensities are already
    % generated.
    STL1_4state_escape = ...
      STL1_4state_escape.formPropensitiesGeneral('STL1_4state_escape');

    % Set the initial populations:
    STL1_4state_escape.initialCondition = [1;0;0;0;0];

    % Set the times at which distributions will be computed:
    STL1_4state_escape.tSpan = linspace(0,100,200);

    % Solve for time for mRNA to reach 100:
    STL1_4state_escape.fspOptions.escapeSinks.f = {'mRNA'};
    STL1_4state_escape.fspOptions.escapeSinks.b = 100;    
    [~,~,STL1_4state_escape] = STL1_4state_escape.solve;

    % Plot the CDF and PDF
    STL1_4state_escape.plotFSP(STL1_4state_escape.Solutions, [],...
        "escapeTimes", [], [], {'linewidth',3}, XLim=[0,50],...
        TitleFontSize=24, Title="4-state STL1 (mRNA)",...
        Colors=[0.23,0.67,0.2], LegendLocation="southeast");
