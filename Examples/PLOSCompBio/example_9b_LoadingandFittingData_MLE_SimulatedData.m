%% example_9b_LoadingandFittingData_MLE_SimulatedData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying Model yeast data 
%   * Maximize the likelihood L(D|theta) and use the maximum likelihood
%     estimate (MLE) to fit the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the two models from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, and fit simulated data from 
% example_8b_LoadingandFittingData_SimulatingData

%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels
% example_4_SolveSSITModels_FSP
% example_8b_LoadingandFittingData_SimulatingData

% View model summaries:
Model_sim_data.summarizeModel
STL1_sim_data.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit simulated data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make copies of the models for MLE:
Model_MLE_sim = Model_sim_data;
STL1_MLE_sim = STL1_sim_data;

% Set fitOptions, with the maximum allowable number of iterations to fit
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
Modelpars_sim = cell2mat(Model_MLE_sim.parameters(1:4,2));
STL1pars_sim = cell2mat(STL1_MLE_sim.parameters(1:7,2));

%% Compute the MLE
[Modelpars_sim,Model_likelihood] = ...
 Model_MLE_sim.maximizeLikelihood(Modelpars_sim,fitOptions);

[STL1pars_sim,STL1_likelihood] = ...
 STL1_MLE_sim.maximizeLikelihood(STL1pars_sim,fitOptions);

% Update Model parameters
for j=1:length(Modelpars_sim)
    Model_MLE_sim.parameters{j,2} = Modelpars_sim(j);
end
for k=1:length(STL1pars_sim)
    STL1_MLE_sim.parameters{k,2} = STL1pars_sim(k);
end

% Make plots of the model parameter fits from the MLEs
Model_MLE_sim.makeFitPlot
STL1_MLE_sim.makeFitPlot