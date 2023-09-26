%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
clear 
close all
clc
addpath('../CommandLine');
addpath('../EricModel/DUSP1_GR_dataframes');

%% Run simulated experiment sampling times
timeMatrix = [0,90,180];
NCells = 100;
dataFileName =  'DUSP1_3hr_Dex_100nM_total.csv'; 
[simData] = sampleExperimentSim(dataFileName,timeMatrix,NCells);