%% example_8_LoadingandFittingData_DataLoading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Data loading and handling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels and 
% example_4_SolveSSITModels_FSP
%clear
%close all
addpath(genpath('../../src'));

%example_1_CreateSSITModels 
%example_4_SolveSSITModels_FSP

% View model summary:
STL1_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load time-varying STL1 yeast data and filter by experimental 
%% conditions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our model:
STL1_data = STL1_FSP;

% Load the experimental data, matching the species name of the model  
% ('mRNA') to the appropriate column in the data file ('rna') and filter 
% it so that the only data loaded into the model is from 'Replica' = 1 and
% 'Condition' = '0.2M_NaCl_Step' (which, in the latter case, is all the 
% data):

STL1_data = STL1_data.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                               {'mRNA','RNA_STL1_total_TS3Full'},...
                               {'Replica',1;'Condition','0.2M_NaCl_Step'});

% This plot is unnecessary, as the model parameters have not been fit to
% the data yet.  However, it illustrates the improvement to come later:
%STL1_data.makeFitPlot