%% example_8_LoadingandFittingData_DataLoading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting Hog1-MAPK pathway STL1 yeast data 
%   * Data loading and handling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels and 
% example_4_SolveSSITModels_FSP
% clear
% close all
addpath(genpath('../src'));

% example_1_CreateSSITModels 
% example_4_SolveSSITModels_FSP

%% Load pre-computed FSP solutions:
% load('example_4_SolveSSITModels_FSP.mat')

% View model summaries:
STL1_4state_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Hog1-MAPK pathway STL1 yeast data and filter by experimental 
%% conditions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our 4-state, time-varying STL1 model:
STL1_4state_data = STL1_4state_FSP;

% Load the experimental data, matching the species name of the model  
% ('mRNA') to the appropriate column in the data file ('rna') and filter 
% it so that the only data loaded into the model is from 'Replica' = 1 and
% 'Condition' = '0.2M_NaCl_Step':
 
STL1_4state_data = ...
    STL1_4state_data.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                             {'mRNA','RNA_STL1_total_TS3Full'},...
                             {'Replica',2;'Condition','0.2M_NaCl_Step'});

STL1_4state.plotFits([], "all", [], {'linewidth',2}, ...
    struct('VarianceType',"STD", 'Title',"STL1 fit"));


% This plot is unnecessary, as the model parameters have not been fit to
% the data yet.  However, it illustrates the improvement to come later:
STL1_4state_data.makeFitPlot

%% Save models with loaded data
saveNames = unique({
    'STL1_4state_data'
    });
    
save('example_8_LoadingandFittingData',saveNames{:})


