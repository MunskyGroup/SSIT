%% example_1b_CreateSSITModels_SimulatingData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries:
%clear
%close all
addpath(genpath('../../'));

%example_1_CreateSSITModels
%example_4_SolveSSITModels_FSP

% Create a copy of the STL1 model for simulation:
STL1_sim = STL1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate simulated time-varying STL1 yeast data for fitting later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of independent data sets to generate:
STL1_sim.ssaOptions.Nexp = 2;  

% Number of cells to include at each time point for each data set:
STL1_sim.ssaOptions.nSimsPerExpt = 200;

% Generate and save data:
dataTable = STL1_sim.sampleDataFromFSP(STL1_FSPsoln,'data/STL1_sim.csv'); 

% Plot data as histograms:
for i = 1:4
    subplot(2,2,i)  % Switch to current subplot
    histogram(dataTable.exp1_s3(dataTable.time==STL1.tSpan(i)),30,...
                                  "DisplayStyle","stairs")
end 