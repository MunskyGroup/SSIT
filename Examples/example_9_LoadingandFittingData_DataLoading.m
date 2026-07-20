%% SSIT/Examples/example_9_LoadingandFittingData_DataLoading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.3.1: Loading and fitting time-varying STL1 yeast data
%   * Data loading and handling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels and 
% example_4_SolveSSITModels_FSP
% clear
% close all

% example_1_CreateSSITModels 
% example_4_SolveSSITModels_FSP

%% Load pre-computed FSP solutions:
load('example_4_SolveSSITModels_FSP.mat')

% View model summaries:
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Hog1-MAPK pathway STL1 yeast data and filter by experimental 
%% conditions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the experimental data, matching the species name of the model  
% ('mRNA') to the appropriate column in the data file 
% ('RNA_STL1_total_TS3Full') and filter it so that the only data loaded 
% into the model is from 'Replica' = 1 and 'Condition' = '0.2M_NaCl_Step':

Model = Model.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                           {'mRNA','RNA_STL1_total_TS3Full'},...
                           {'Replica',1;'Condition','0.2M_NaCl_Step'});

STL1 = STL1.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                           {'mRNA','RNA_STL1_total_TS3Full'},...
                           {'Replica',1;'Condition','0.2M_NaCl_Step'});
 
STL1_4state = ...
    STL1_4state.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                           {'mRNA','RNA_STL1_total_TS3Full'},...
                           {'Replica',1;'Condition','0.2M_NaCl_Step'});

% These plots are unnecessary, as the model parameters have not been fit
% to the data yet.  However, it illustrates the improvement to come later:

Model.plotFits(plotType="all", lineProps={'linewidth',2},...
    Title='Bursting Gene', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

STL1.plotFits(plotType="all", lineProps={'linewidth',2},...
    Title='STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

STL1_4state.plotFits(plotType="all", lineProps={'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

%% Save models with loaded data
saveNames = unique({
    'Model'
    'STL1'
    'STL1_4state'
    });
    
save('example_9_LoadingandFittingData_DataLoading',saveNames{:}) 