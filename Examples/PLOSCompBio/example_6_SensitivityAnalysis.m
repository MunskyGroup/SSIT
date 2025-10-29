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

% Make a copy of the bursting gene model solved by FSP for sensitivity 
% analysis:
Model_sens = Model_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
Model_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem:
[Model_sensSoln,Model_bounds] = Model_sens.solve(Model_FSPsoln.stateSpace);

% Plot the results from the sensitivity analysis:
fig1 = figure(1);clf; set(fig1,'Name','Marginal Sensitivity, offGene');
fig2 = figure(2);clf; set(fig2,'Name','Marginal Sensitivity, onGene');
fig3 = figure(3);clf; set(fig3,'Name','Marginal Sensitivity, mRNA');
Model_sens.makePlot(Model_sensSoln,'marginals',[],false,...
                    [fig1,fig2,fig3],{'b','linewidth',2})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve sensitivities of the time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the time-varying STL1 yeast model solved by FSP for 
% sensitivity analysis:
STL1_sens = STL1_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
STL1_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem: 
[STL1_sensSoln,STL1_bounds] = STL1_sens.solve(STL1_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis:
fig4 = figure(4);clf; set(fig4,'Name','Marginal Sensitivity, offGene');
fig5 = figure(5);clf; set(fig5,'Name','Marginal Sensitivity, onGene');
fig6 = figure(6);clf; set(fig6,'Name','Marginal Sensitivity, mRNA');
STL1_sens.makePlot(STL1_sensSoln,'marginals',[],false,...
                   [fig4,fig5,fig6],{'b','linewidth',2})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve sensitivities of the 4-state time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the time-varying STL1 yeast model solved by FSP for 
% sensitivity analysis:
STL1_4state_sens = STL1_4state_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
STL1_4state_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem: 
[STL1_4state_sensSoln,STL1_4state_bounds] = ...
    STL1_4state_sens.solve(STL1_4state_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis:
% fig7 = figure(7);clf; set(fig7,'Name','Marginal Sensitivity, s1');
% fig8 = figure(8);clf; set(fig8,'Name','Marginal Sensitivity, s2');
% fig9 = figure(9);clf; set(fig9,'Name','Marginal Sensitivity, s3');
% fig10 = figure(10);clf; set(fig10,'Name','Marginal Sensitivity, s4');
% fig11 = figure(11);clf; set(fig11,'Name','Marginal Sensitivity, mRNA');
% STL1_4state_sens.makePlot(STL1_4state_sensSoln,'marginals',[],false,...
%                    [fig7,fig8,fig9,fig10,fig11],{'b','linewidth',2})

%% Plot the 
STL1_4state_sens.plotFSP(STL1_4state_sensSoln, STL1_4state_sens.species,...
        'fspSens', [], [], {'linewidth',3}, AxisLabelSize=18,...
        TickLabelSize=18)

%% Save models & sensitivities
saveNames = unique({'Model_sens'
    'Model_sensSoln'
    'Model_bounds'
    'STL1_sens'
    'STL1_sensSoln'
    'STL1_bounds'
    'STL1_4state_sens'
    'STL1_4state_sensSoln'
    'STL1_4state_bounds'
    });

    
    %% Specify a boundary for the escape calculation
    % Calculate the time until the mRNA concentration reaches 5
    Model_escape.fspOptions.escapeSinks.f = {'mRNA'};
    Model_escape.fspOptions.verbose = false;
    Model_escape.fspOptions.escapeSinks.b = 5;
    Model_escape = Model_escape.formPropensitiesGeneral('Model_escape');
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
    STL1_escape = STL1_escape.formPropensitiesGeneral('STL1_escape');
    
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
    STL1_4state_escape = ...
        STL1_4state_escape.formPropensitiesGeneral('STL1_4state_escape');
    
    % Solve for the escape time:
    STL1_4state_escape.fspOptions.escapeSinks.f = {'g4'}
    STL1_4state_escape.fspOptions.escapeSinks.b = [0.5];
    [STL1_4state_fspSoln_escape,STL1_4state_escape.fspOptions.bounds] = ...
        STL1_4state_escape.solve;

    % Plot the CDF and PDF
    fig3 = figure(3); clf; set(fig3,'Name','4-state STL1');
    STL1_4state_escape.makePlot(STL1_4state_fspSoln_escape,...
                                'escapeTimes',[],false,fig3)
    legend(STL1_4state_escape.fspOptions.escapeSinks.f)
