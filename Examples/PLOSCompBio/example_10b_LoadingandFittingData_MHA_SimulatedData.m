%% example_10b_LoadingandFittingData_MHA_SimulatedData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying Model yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the two models from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, simulated data from 
% example_1b_CreateSSITModels_SimulatingData, data loaded in 
% example_8b_LoadingandFittingData_SimulatedDataLoading, and MLE computed 
% in example_9b_LoadingandFittingData_MLE_SimulatedData
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_1b_CreateSSITModels_SimulatingData
% example_8b_LoadingandFittingData_SimulatedDataLoading
% example_9b_LoadingandFittingData_MLE_SimulatedData

% View model summaries:
Model_MLE_sim.summarizeModel
STL1_MLE_sim.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Metropolis-Hastings to sample uncertainty 
%   (and improve model parameter fit to data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of our models for Metropolis-Hastings (MH):
Model_MH_sim = Model_MLE_sim;
STL1_MH_sim = STL1_MLE_sim;

    % Adjust proposal width scale (the default proposal distribution in 
    % SSIT is "@(x)x+0.1*randn(size(x))", which leads to low acceptance in
    % this case. A rule of thumb, depending on the data being analyzed, is 
    % to aim for an MH acceptance ratio of around 0.3-0.4 (some say 25%):
    proposalWidthScale = 0.005;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

% Set MH runtime options - number of samples, burnin, thin sample set,
% etc.
MHOptions.numberOfSamples = 1000;
MHOptions.burnin = 100;
MHOptions.thin = 1;

% Run Metropolis-Hastings: 
[Model_MH_sim_pars,~,Model_sim_MHResults] = ...
   Model_MH_sim.maximizeLikelihood([], MHOptions, 'MetropolisHastings');
[STL1_MH_sim_pars,~,STL1_sim_MHResults] = ...
   STL1_MH_sim.maximizeLikelihood([], MHOptions, 'MetropolisHastings');

% Store MH parameters in model:
Model_MH_sim.parameters([1:4],2) = num2cell(Model_MH_sim_pars);
STL1_MH_sim.parameters([1:7],2) = num2cell(STL1_MH_sim_pars);

% Plot results:
Model_MH_sim.plotMHResults(Model_sim_MHResults,[],'log',[])
Model_MH_sim.makeFitPlot
STL1_MH_sim.plotMHResults(STL1_sim_MHResults,[],'log',[])
STL1_MH_sim.makeFitPlot

%% Specify Bayesian Prior and fit
% Specify the prior as log-normal distribution with wide uncertainty
mu_log10 = [-1,-1,1,1];  % Prior log-mean
sig_log10 = 2*ones(1,4);  % Prior log-standard deviation
Model_MH_sim.fittingOptions.logPrior = @(x)-sum((log10(x)-mu_log10).^2./...
                                        (2*sig_log10.^2));

% Choose parameters to search (in this case, all 7 model parameters):
Model_MH_sim.fittingOptions.modelVarsToFit = [1:4];

% Create initial parameter guess:
Model_MH_sim_pars = [Model_MH_sim.parameters{:,2}]; 

 % Fit to maximize likelihood:
Model_MH_sim_pars = Model_MH_sim.maximizeLikelihood(Model_MH_sim_pars);

% Update new parameters:
Model_MH_sim.parameters(:,2) = num2cell(Model_MH_sim_pars);

% Plot fitting results
Model_MH_sim.makeFitPlot  

%% STL1 %%
mu_log10 = [-1,1,1,1,1,-2,-2];    % Prior log-mean
sig_log10 = 2*ones(1,7);          % Prior log-standard deviation
STL1_MH_sim.fittingOptions.logPrior = @(x)-sum((log10(x)-mu_log10).^2./...
                                       (2*sig_log10.^2));

% Choose parameters to search (in this case, all 7 model parameters):
STL1_MH_sim.fittingOptions.modelVarsToFit = [1:7];

% Create initial parameter guess:
STL1_MH_sim_pars = [STL1_MH_sim.parameters{:,2}]; 

 % Fit to maximize likelihood:
STL1_MH_sim_pars = STL1_MH_sim.maximizeLikelihood(STL1_MH_sim_pars);

% Update new parameters:
STL1_MH_sim.parameters(:,2) = num2cell(STL1_MH_sim_pars);

% Plot fitting results
STL1_MH_sim.makeFitPlot  