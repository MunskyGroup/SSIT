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
    proposalWidthScale = 0.01;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

% Set MH runtime options - number of samples, burnin, thin sample set,
% etc.
MHOptions.numberOfSamples = 3000;
MHOptions.burnin = 300;
MHOptions.thin = 3;

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
mu_log10 = [-1,1,-1,1,1,-2,-2];   % Prior log-mean
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

%% Iterating between MLE and MH
% Let's run a few rounds of MLE and MH to see if we can get better
% convergence.
Model_MH_sim.parameters(:,2) = num2cell(Model_MH_sim_pars);
for i=1:3
    % Maximize likelihood
    Model_MH_sim_pars = Model_MH_sim.maximizeLikelihood([],fitOptions);    
    % Update parameters in the model:
    Model_MH_sim.parameters(:,2) = num2cell(Model_MH_sim_pars);
    
    % Compute FIM
    % Set solutions scheme to FSP Sensitivity
    Model_MH_sim.solutionScheme = 'fspSens'; 
    % Solve the sensitivity problem
    [Model_sensSoln] = Model_MH_sim.solve;  
    Model_fimResults = Model_MH_sim.computeFIM(Model_sensSoln.sens,'log');
    Model_FIMlog = Model_MH_sim.evaluateExperiment(Model_fimResults,...
                                              Model_MH_sim.dataSet.nCells);

    % Run Metropolis-Hastings    
    % Adjusted proposal dist. covariance.
    Model_covLogMod = (Model_FIMlog{1} + ...
                       diag(size(Model_FIMlog{1},1)))^(-1); 
    proposalWidthScale = 1;
    Model_MH_simOptions.proposalDistribution = ...
     @(x)mvnrnd(x,proposalWidthScale*(Model_covLogMod+Model_covLogMod')/2);
    [Model_MH_sim_pars,Model_likelihood_FIM,Model_chainResults_FIM] = ...
        Model_MH_sim.maximizeLikelihood([],Model_MH_simOptions,...
                                         'MetropolisHastings');
    % Update parameters in the model:
    Model_MH_sim.parameters(:,2) = num2cell(Model_MH_sim_pars);
end
Model_MH_sim.plotMHResults(Model_chainResults_FIM,Model_FIMlog);
Model_MH_sim.makeFitPlot
%%
STL1_MH_sim.parameters(:,2) = num2cell(STL1_MH_sim_pars);
for i=1:3
    % Maximize likelihood
    STL1_MH_sim_pars = STL1_MH_sim.maximizeLikelihood([],fitOptions);    
    % Update parameters in the model:
    STL1_MH_sim.parameters(:,2) = num2cell(STL1_MH_sim_pars);
    
    % Compute FIM
    % Set solutions scheme to FSP Sensitivity
    STL1_MH_sim.solutionScheme = 'fspSens'; 
    % Solve the sensitivity problem
    [STL1_sensSoln] = STL1_MH_sim.solve;  
    STL1_fimResults = STL1_MH_sim.computeFIM(STL1_sensSoln.sens,'log');
    STL1_FIMlog = STL1_MH_sim.evaluateExperiment(STL1_fimResults,...
                                               STL1_MH_sim.dataSet.nCells);

    % Run Metropolis-Hastings    
    % Adjusted proposal dist. covariance.
    STL1_covLogMod = (STL1_FIMlog{1} + diag(size(STL1_FIMlog{1},1)))^(-1); 
    proposalWidthScale = 0.001;
    STL1_MH_simOptions.proposalDistribution = ...
       @(x)mvnrnd(x,proposalWidthScale*(STL1_covLogMod+STL1_covLogMod')/2);
    [STL1_MH_sim_pars,STL1_likelihood_FIM,STL1_chainResults_FIM] = ...
        STL1_MH_sim.maximizeLikelihood([],STL1_MH_simOptions,...
                                        'MetropolisHastings');
    % Update parameters in the model:
    STL1_MH_sim.parameters(:,2) = num2cell(STL1_MH_sim_pars);
end
STL1_MH_sim.plotMHResults(STL1_chainResults_FIM,STL1_FIMlog);
STL1_MH_sim.makeFitPlot


%% Evaluating the MH results
% Here we will generate three plots.  The first one will show the
% Likelihood function as we search over parameter space.  In order to get a
% good estimate of the parameter uncertainty, we want this to quickly reach
% ts maximum value and then to fluctuate around that value for a
% significant amount of time.  If you see that it is still incereasing, you
% know that the fit has not yet converged.
figure
plot(Model_chainResults_FIM.mhValue)
title('MH Convergence')
xlabel('Iteration Number')
ylabel('LogLikelihood')

figure
plot(STL1_chainResults_FIM.mhValue)
title('MH Convergence')
xlabel('Iteration Number')
ylabel('LogLikelihood')

%% Compute FIM
%% Model:
% Set solutions scheme to FSP Sensitivity
Model_MH_sim.solutionScheme = 'fspSens'; 
% Solve the sensitivity problem
[Model_sensSoln] = Model_MH_sim.solve;  
Model_fimResultsLog = Model_MH_sim.computeFIM(Model_sensSoln.sens,'log');
Model_FIMlog = Model_MH_sim.evaluateExperiment(Model_fimResultsLog,Model_MH_sim.dataSet.nCells);
Model_fimResults= Model_MH_sim.computeFIM(Model_sensSoln.sens,'lin');
Model_FIMlin = Model_MH_sim.evaluateExperiment(Model_fimResults,Model_MH_sim.dataSet.nCells);

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.

% Choose which parameters to compare.
Q = [3,4];  
subplot(1,3,2)

% Plot uncertainty in linear scale
Model_MH_sim_FIM.makeMleFimPlot(exp(Model_chainResults_FIM.mhSamples)',...
                                Model_FIMlin{1},Q,0.95,1); hold on
title('Posterior -- Linear Scale')
xlabel(Model_MH_sim_FIM.parameters{Q(1),1})
ylabel(Model_MH_sim_FIM.parameters{Q(2),1})

% Plot uncertainty in log scale
subplot(1,3,3)
Model_MH_sim_FIM.makeMleFimPlot(Model_chainResults_FIM.mhSamples',...
                                Model_FIMlog{1},Q,0.95,1); hold on
title('Posterior -- Natural Log Scale')
xlabel(['log ',Model_MH_sim_FIM.parameters{Q(1),1}])
ylabel(['log ',Model_MH_sim_FIM.parameters{Q(2),1}]) 

%% STL1:
% Set solutions scheme to FSP Sensitivity
STL1_MH_sim.solutionScheme = 'fspSens'; 
% Solve the sensitivity problem
[STL1_sensSoln] = STL1_MH_sim.solve;  
STL1_fimResultsLog = STL1_MH_sim.computeFIM(STL1_sensSoln.sens,'log');
STL1_FIMlog = STL1_MH_sim.evaluateExperiment(STL1_fimResultsLog,STL1_MH_sim.dataSet.nCells);
STL1_fimResults= STL1_MH_sim.computeFIM(STL1_sensSoln.sens,'lin');
STL1_FIMlin = STL1_MH_sim.evaluateExperiment(STL1_fimResults,STL1_MH_sim.dataSet.nCells);

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.

% Choose which parameters to compare.
Q = [3,4];  
subplot(1,3,2)

% Plot uncertainty in linear scale
STL1_MH_sim_FIM.makeMleFimPlot(exp(STL1_chainResults_FIM.mhSamples)',...
                                STL1_FIMlin{1},Q,0.95,1); hold on
title('Posterior -- Linear Scale')
xlabel(STL1_MH_sim_FIM.parameters{Q(1),1})
ylabel(STL1_MH_sim_FIM.parameters{Q(2),1})

% Plot uncertainty in log scale
subplot(1,3,3)
STL1_MH_sim_FIM.makeMleFimPlot(STL1_chainResults_FIM.mhSamples',...
                                STL1_FIMlog{1},Q,0.95,1); hold on
title('Posterior -- Natural Log Scale')
xlabel(['log ',STL1_MH_sim_FIM.parameters{Q(1),1}])
ylabel(['log ',STL1_MH_sim_FIM.parameters{Q(2),1}]) 

%% Effective Sample Size
% For the MH, it is important to get a sense of how well it has sampled the
% posterior.  For this, we determine the effective sample size (i.e., the
% number of effectively indpendent samples within the MH chain).  This is
% found by examining the autocorrelation of the parameter chain to figure
% out the number of steps needed for correlations to decay and then divide
% the total number of steps by the de-correlation step.
figure
ipar = 4;
ac = xcorr(STL1_chainResults_FIM.mhSamples(:,ipar)-mean(STL1_chainResults_FIM.mhSamples(:,ipar)),'normalized');
ac = ac(size(STL1_chainResults_FIM.mhSamples,1):end);
plot(ac,'LineWidth',3)
N = size(STL1_chainResults_FIM.mhSamples,1);
tau = 1+2*sum(abs(ac(2:N/5)));
Neff = N/tau