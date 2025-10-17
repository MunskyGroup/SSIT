%% example_5_SolveSSITModels_EscapeTimes 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Solve a first-passage time problem (escape times)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels
% clear
% close all
% addpath(genpath('../../'));

% example_1_CreateSSITModels

% Load the models created in example_1_CreateSSITModels
% load('example_1_CreateSSITModels.mat')

% View model summaries:
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

% Set the times at which distributions will be computed:
Model.tSpan = linspace(0,20,200);
STL1.tSpan = linspace(0,20,200);
STL1_4state.tSpan = linspace(0,50,200);

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
    Model_escape.formPropensitiesGeneral('Model_escape');
    [fspSoln_escape_1,Model_escape.fspOptions.bounds] = Model_escape.solve;

    % Plot the CDF and PDF
    fig1 = figure(1); clf; set(fig1,'Name','Bursting Gene');
    Model_escape.makePlot(fspSoln_escape_1,'escapeTimes',[],false,fig1)
    legend(Model_escape.fspOptions.escapeSinks.f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve escape times for the time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1:
    % Create a copy of the time-varying STL1 yeast model:
    STL1_escape = STL1;
    STL1_escape.formPropensitiesGeneral('STL1_escape');
    
    % Solve for the escape time:
    STL1_escape.fspOptions.escapeSinks.f = {'offGene';'onGene'};
    STL1_escape.fspOptions.escapeSinks.b = [0.5;0.5];
    [STL1_fspSoln_escape,STL1_escape.fspOptions.bounds] = ...
        STL1_escape.solve;

    % Plot the CDF and PDF
    fig2 = figure(2); clf; set(fig2,'Name','STL1');
    STL1_escape.makePlot(STL1_fspSoln_escape,'escapeTimes',[],false,fig2)
    legend(STL1_escape.fspOptions.escapeSinks.f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve escape times for the 4-state time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4-state STL1:
    % Create a copy of the time-varying STL1 yeast model:
    STL1_4state_escape = STL1_4state;
    STL1_4state_escape.formPropensitiesGeneral('STL1_4state_escape');
    
    % Solve for the escape time:
    STL1_4state_escape.fspOptions.escapeSinks.f = {'g2+g3'}
    STL1_4state_escape.fspOptions.escapeSinks.b = [0.5];
    [STL1_4state_fspSoln_escape,STL1_4state_escape.fspOptions.bounds] = ...
        STL1_4state_escape.solve;

    % Plot the CDF and PDF
    STL1_4state_escape.plotFSP(STL1_4state_fspSoln_escape, [],...
        "escapeTimes", 1:5:200, [], {'linewidth',2},...
        Title="First-Passage Statistics", XLabel="Time", XLim=[0 50]);