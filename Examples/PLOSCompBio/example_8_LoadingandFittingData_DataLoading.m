%% example_8_LoadingandFittingData_DataLoading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Data loading and handling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels 
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  

% View model summary:
STL1_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load time-varying STL1 yeast data and filter by experimental 
%% condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our model:
STL1_data = STL1_FSP;

%% Load the experimental data, matching the species name of the model  
%% ('mRNA') to the appropriate column in the data file ('rna') and filter 
%% it so that only data from experiment 1, replica 1 is loaded:

% STL1_data = STL1_data.loadData('data/Hog1_exp_rep_CY5_combined.csv',...
%                               {'mRNA','rna'},{'experiment',1;'replica',1});
STL1_data = STL1_data.loadData('data/STL1.csv',...
                               {'mRNA','rna'});

% This plot is unnecessary, as the model parameters have not been fit to
% the data yet.  However, it illustrates the improvement to come later:
STL1_data.makeFitPlot