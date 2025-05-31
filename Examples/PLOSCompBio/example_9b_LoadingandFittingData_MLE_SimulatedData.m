%% example_9b_LoadingandFittingData_MLE_SimulatedData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Maximize the likelihood L(D|theta) and use the maximum likelihood
%     estimate (MLE) to fit the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the simulated STL1 data from example_1b_CreateSSITModels and 
% example_8b_LoadingandFittingData_SimulatedDataLoading
%clear
%close all
addpath(genpath('../../src'));

% example_1b_CreateSSITModels_SimulatingData 
% example_8b_LoadingandFittingData_SimulatedDataLoading

% View model summary:
STL1_sim_data.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit simulated data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of STL1_data for MLE:
STL1_MLE_sim = STL1_sim_data;

% Set fitOptions, with the maximum allowable number of iterations to fit
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
STL1pars_sim = cell2mat(STL1_MLE_sim.parameters(1:7,2));

%% Compute the MLE
[STL1pars_sim,STL1_likelihood] = ...
 STL1_MLE_sim.maximizeLikelihood(STL1pars_sim,fitOptions);

% Update STL1 parameters
for j=1:length(STL1pars_sim)
    STL1_MLE_sim.parameters{j,2} = STL1pars_sim(j);
end

% Make plots of the model parameter fits from the MLEs
STL1_MLE_sim.makeFitPlot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's tinker with the starting parameters of the STL1 Model and try 
%% again:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of STL1_MLE_sim for refitting:
STL1_MLE_sim_refit = STL1_MLE_sim;

STL1_MLE_sim_refit.parameters = ({'koff',30;'kr',100;'gr',0.005; ...
                              'a0',0.01;'a1',1;'r1',0.4;'r2',0.1});

STL1_MLE_sim_refit.tSpan = 0:5:60;

% Have SSIT approximate the steady state for initial distribution
STL1_MLE_sim_refit.fspOptions.initApproxSS =true;
STL1_MLE_sim_refit = ...
    STL1_MLE_sim_refit.formPropensitiesGeneral('STL1_MLE_sim_refit',true);

% Solve the FSP analysis
[STL1_MLE_sim_refit_FSPsoln,STL1_MLE_sim_refit.fspOptions.bounds] = ...
    STL1_MLE_sim_refit.solve;  

% Format new parameters
STL1pars_sim_refit = cell2mat(STL1_MLE_sim_refit.parameters(1:7,2));

% Maximize the likelihood
[STL1pars_sim_refit,STL1_likelihood_refit] = ...
    STL1_MLE_sim_refit.maximizeLikelihood(STL1pars_sim_refit, fitOptions);

% Update STL1Real parameters
for j=1:length(STL1pars_sim_refit)
    STL1_MLE_sim_refit.parameters{j,2} = STL1pars_sim_refit(j);
end

% Make plots of the new model parameter fits from the MLEs
STL1_MLE_sim_refit.makeFitPlot