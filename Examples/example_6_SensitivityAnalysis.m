%% SSIT/Examples/example_6_SensitivityAnalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Sensitivity analysis 
%   * Compute model sensitivities to small changes in the model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels and  
% computed FSP solutions from example_4_SolveSSITModels_FSP

% clear
% close all

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
[~,~,Model_sens] = Model_sens.solve(Model_FSP.Solutions.stateSpace);

% Plot the results from the sensitivity analysis:
Model_sens.plotFSP(Model_sens.Solutions,Model_FSP.species(3),'sens',40,...
   [], {'linewidth',3}, AxisLabelSize=12, TickLabelSize=12, XLim=[0,10],...
   TitleFontSize=22, Title="Bursting Gene (mRNA)", Colors=[0.93,0.69,0.13])

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
[~,~,STL1_sens] = STL1_sens.solve(STL1_FSP.Solutions.stateSpace); 

% Plot the results from the sensitivity analysis:
STL1_sens.plotFSP(STL1_sens.Solutions,STL1_FSP.species(3),'sens',40,[],...
    {'linewidth',3}, AxisLabelSize=12, TickLabelSize=12, XLim=[0,10],...
    Title="STL1 (mRNA)", TitleFontSize=22, Colors=[0.93,0.69,0.13])

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
[~,~,STL1_4state_sens] = STL1_4state_sens.solve;

% Plot the results from the sensitivity analysis
STL1_4state_sens.plotFSP(STL1_4state_sens.Solutions,...
    STL1_4state_sens.species(5), 'sens', 40, [],...
    {'linewidth',3}, Colors=[0.23,0.67,0.2], AxisLabelSize=15,...
    TickLabelSize=12, XLim=[0,100],...
    Title="4-state STL1 (t=25)", TitleFontSize=24)


%% Save models & sensitivities
saveNames = unique({ ...
    'Model_sens'
    'STL1_sens'
    'STL1_4state_sens'
    });
    
save('example_6_SensitivityAnalysis',saveNames{:})