%% example_8b_LoadingandFittingData_SimulatedDataLoading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting simulated time-varying STL1 yeast data 
%   * Data loading and handling 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the time-varying STL1 yeast model from example_1_CreateSSITModels, 
% FSP solutions from example_4_SolveSSITModels_FSP, and 
% simulated STL1 data from example_1b_CreateSSITModels_SimulatingData 
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels
% example_4_SolveSSITModels_FSP
% example_1b_CreateSSITModels_SimulatingData

% View model summary:
Model_FSP.summarizeModel
STL1_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load time-varying STL1 yeast data and filter by experimental 
%% conditions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make new copies of our model:
Model_sim_data = Model_FSP;
STL1_sim_data = STL1_FSP;

% Load the simulated data, matching the species name of the model  
% (e.g., the model species 'offGene' is column 'exp1_s1' in the file):

Model_sim_data = Model_sim_data.loadData('data/STL1_sim.csv',...
                {'offGene','exp1_s1';'onGene','exp1_s2';'mRNA','exp1_s3'});

STL1_sim_data = STL1_sim_data.loadData('data/STL1_sim.csv',...
                {'offGene','exp1_s1';'onGene','exp1_s2';'mRNA','exp1_s3'});

% This plot is unnecessary, as the model parameters have not been fit to
% the data yet.  However, it illustrates the improvement to come later:
Model_sim_data.makeFitPlot
STL1_sim_data.makeFitPlot