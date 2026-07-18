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
load('example_4_SolveSSITModels_FSP.mat')

% View model summaries
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Solve sensitivities of the bursting gene model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve FSP sensitivities
% Solve the sensitivity problem:
Model = Model.solve(solver='fspSens');

% Plot the results from the sensitivity analysis:
Model.plotFSP(speciesNames=Model.species(3), plotType='sens',...
   indTimes=40, lineProps={'linewidth',3}, AxisLabelSize=12,...
   TickLabelSize=12, XLim=[0,10], TitleFontSize=22,...
   Title="Bursting Gene", Colors=[0.93,0.69,0.13])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve sensitivities of the time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve FSP sensitivities
% Solve the sensitivity problem: 
STL1 = STL1.solve(solver='fspSens'); 

% Plot the results from the sensitivity analysis:
STL1.plotFSP(speciesNames=STL1.species(3), plotType='sens',...
    indTimes=40, lineProps={'linewidth',3}, AxisLabelSize=12,...
    TickLabelSize=12, XLim=[0,10], Title="STL1",...
    TitleFontSize=22, Colors=[0.93,0.69,0.13])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve sensitivities of the 4-state time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve FSP sensitivities
% Solve the sensitivity problem:
STL1_4state = STL1_4state.solve(solver='fspSens');

% Plot the results from the sensitivity analysis
STL1_4state.plotFSP(speciesNames=STL1_4state.species(5),...
    plotType='sens', indTimes=40, lineProps={'linewidth',3},...
    Colors=[0.23,0.67,0.2], AxisLabelSize=15, TickLabelSize=12,...
    XLim=[0,100], Title="4-state STL1 (t=25)", TitleFontSize=24)


%% Save models & sensitivities
saveNames = unique({ ...
    'Model'
    'STL1'
    'STL1_4state'
    });
    
save('example_6_SensitivityAnalysis',saveNames{:})