%% example_6_SensitivityAnalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.3: Sensitivity analysis and Fisher Information Matrix (Part I)
%   * Compute model sensitivities to small changes in the model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels and  
% computed FSP solutions from example_4_SolveSSITModels_FSP

% clear
% close all
% addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP

%% Load pre-computed FSP solutions:
% load('example_4_SolveSSITModels_FSP.mat')

% View model summaries
Model_FSP.summarizeModel
STL1_FSP.summarizeModel
STL1_4state_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Solve sensitivities of the bursting gene model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the bursting gene model solved by FSP for sensitivity 
% analysis:
Model_sens = Model_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
Model_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem:
[Model_sensSoln,Model_bounds] = Model_sens.solve(Model_FSPsoln.stateSpace);

% Plot the results from the sensitivity analysis:
Model_sens.plotFSP(Model_sensSoln, Model_FSP.species(3), 'sens', 20, [],...
    {'linewidth',3}, AxisLabelSize=12, TickLabelSize=12)

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
STL1_sens.plotFSP(STL1_sensSoln, STL1_FSP.species(3), 'sens', 20, [],...
    {'linewidth',3}, AxisLabelSize=12, TickLabelSize=12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve sensitivities of the 4-state time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the 4-state time-varying STL1 yeast model solved by FSP:
STL1_4state_sens = STL1_4state_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
STL1_4state_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem: 
[STL1_4state_sensSoln,STL1_4state_bounds] = ...
    STL1_4state_sens.solve(STL1_4state_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis:
STL1_4state_sens.plotFSP(STL1_4state_sensSoln,...
    STL1_4state_FSP.species(5), 'sens', 20, [], {'linewidth',3}, ...
    Colors=[0.23,0.67,0.2], AxisLabelSize=12, TickLabelSize=12, ...
    XLim=[0,100])


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
    
save('example_6_SensitivityAnalysis',saveNames{:})