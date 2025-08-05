%% example_10_LoadingandFittingData_MHA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, sensitivities 
% from example_6_SensitivityAnalysis, the FIM from example_7_FIM, data   
% loaded in example_8_LoadingandFittingData_DataLoading, and MLE computed 
% in example_9_LoadingandFittingData_MLE
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_6_SensitivityAnalysis
% example_7_FIM
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE

% View summary of 4-state STL1 model:
STL1_MLE_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Metropolis-Hastings to sample uncertainty 
%   (and improve model parameter fit to data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of our 4-state STL1 model for Metropolis-Hastings (MH):
STL1_4state_MH = STL1_MLE_4state;

    % Adjust proposal width scale (the default proposal distribution in 
    % SSIT is "@(x)x+0.1*randn(size(x))", which leads to low acceptance in
    % this case. A rule of thumb, depending on the data being analyzed, is 
    % to aim for an MH acceptance ratio of around 0.3-0.4 (some say 25%):
    proposalWidthScale = 0.01;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

% Set MH runtime options (number of samples, burnin, thin, etc.):
MHOptions.numberOfSamples = 2000;
MHOptions.burnin = 200;
MHOptions.thin = 2;

% Run Metropolis-Hastings: 
[STL1_MH_pars_4state,~,STL1_MHResults_4state] = ...
    STL1_4state_MH.maximizeLikelihood([], MHOptions, 'MetropolisHastings');

% Store MH parameters in model:
STL1_4state_MH.parameters([1:15],2) = num2cell(STL1_MH_pars_4state);

% Plot results:
STL1_4state_MH.plotMHResults(STL1_MHResults_4state,[],'log',[])
STL1_4state_MH.makeFitPlot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use FIM for Metropolis-Hastings proposal distribution 
%   This sometimes provides faster mixing (convergence), although in our 
%   simple example, the default proposal distribution (above) is fine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our 4-state STL1 model:
STL1_4state_MH_FIM = STL1_MLE_4state;

%% Compute FIM, Run Metropolis Hastings
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [-1,-1,-0.5,-3,-4,-7,-10,-0.5,1,-1,-0.5,1,-5,1,-15];  

% Prior log-standard deviation:
sig_log10 = 2*ones(1,15);      

% Prior:
STL1_4state_MH_FIM.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose parameters to search:
STL1_4state_MH_FIM.fittingOptions.modelVarsToFit = [1:15];

% Create first parameter guess:
STL1_4state_MH_FIM_pars = [STL1_4state_MH_FIM.parameters{:,2}];         

% Compute individual FIMs:
fimResults = STL1_4state_MH_FIM.computeFIM([],'log'); 

% Compute total FIM including effect of prior:
fimTotal = STL1_4state_MH_FIM.evaluateExperiment(fimResults,...
           STL1_4state_MH_FIM.dataSet.nCells,diag(sig_log10.^2)); 

% Choose parameters to search:
STL1_4state_MH_FIM.fittingOptions.modelVarsToFit = [1:15]; 

% Select FIM for free parameters:
FIMfree = fimTotal{1}([1:15],[1:15]); 

% Estimate the covariance using CRLB:
COVfree = (1/2*(FIMfree + FIMfree'))^(-1);  

% Define Metropolis-Hasting settings:
STL1_4state_MH_FIM.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10([1:15])).^2./(2*sig_log10([1:15]).^2));
proposalWidthScale = 0.01;
STL1_4state_MH_FIM_FIMOptions = ...
 struct('proposalDistribution',@(x)mvnrnd(x,proposalWidthScale*COVfree),...
 'numberOfSamples',2000,'burnin',200,'thin',2);

% Run Metropolis Hastings
[STL1_4state_MH_FIM_pars,~,STL1_FIM_MHResults] = ...
    STL1_4state_MH_FIM.maximizeLikelihood([], ...
    STL1_4state_MH_FIM_FIMOptions, 'MetropolisHastings'); 

% Store sampled parameters:
STL1_4state_MH_FIM.parameters([1:15],2) = ...
    num2cell(STL1_4state_MH_FIM_pars);

% Plot MH samples, FIM:
STL1_4state_MH_FIM.plotMHResults(STL1_FIM_MHResults,FIMfree,'log',[])
STL1_4state_MH_FIM.makeFitPlot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use Bayesian priors and iterate between computing MLE and MH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our 4-state STL1 model:
STL1_4state_MH_it = STL1_4state_MH_FIM;

%% Specify Bayesian Prior and fit
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [-1,-1,-0.5,-3,-4,-7,-10,-0.5,1,-1,-0.5,1,-5,1,-15]; 

% Prior log-standard deviation:
sig_log10 = 2*ones(1,15);  

% Prior:
STL1_4state_MH_it.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose parameters to search:
STL1_4state_MH_it.fittingOptions.modelVarsToFit = [1:15]; 

% Create first parameter guess:
STL1_4state_MH_it_pars = [STL1_4state_MH_it.parameters{:,2}];      

% Fit to maximize likelihood:
STL1_4state_MH_it_pars = ...
    STL1_4state_MH_it.maximizeLikelihood(STL1_4state_MH_it_pars); 

% Update new parameters:
STL1_4state_MH_it.parameters(:,2) = num2cell(STL1_4state_MH_it_pars); 

% Plot fitting results:
STL1_4state_MH_it.makeFitPlot  
% You may need to re-run this multiple times until converged.
% I got a MLE of -52,454.1 after a few runs. 


%% Iterating between MLE and MH
%  Running a few rounds of MLE and MH together may improve convergence.

STL1_4state_MH_it.parameters(:,2) = num2cell(STL1_4state_MH_it_pars);
for i=1:3
    % Maximize likelihood:
    STL1_4state_MH_it_pars = STL1_4state_MH_it.maximizeLikelihood([]);    
    % Update parameters in the model:
    STL1_4state_MH_it.parameters(:,2) = num2cell(STL1_4state_MH_it_pars);

    % Run Metropolis-Hastings    
    proposalWidthScale = 0.01;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

    % Set MH runtime options (number of samples, burnin, thin, etc.):
    MHOptions.numberOfSamples = 1000;
    MHOptions.burnin = 100;
    MHOptions.thin = 1;

    % Run Metropolis-Hastings: 
    [STL1_4state_MH_it_pars,~,STL1_4state_MH_it_MHResults] = ...
        STL1_4state_MH_it.maximizeLikelihood([], MHOptions,...
        'MetropolisHastings');
    
    % Store MH parameters in model:
    STL1_4state_MH_it.parameters([1:15],2) = ...
        num2cell(STL1_4state_MH_it_pars);
end
STL1_4state_MH_it.plotMHResults(STL1_4state_MH_it_MHResults);
STL1_4state_MH_it.makeFitPlot

%% Evaluating the MH results
