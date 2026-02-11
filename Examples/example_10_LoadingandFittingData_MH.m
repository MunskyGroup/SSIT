%% SSIT/Examples/example_10_LoadingandFittingData_MH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, data loaded in 
% example_8_LoadingandFittingData_DataLoading, and MLE computed in
% example_9_LoadingandFittingData_MLE
% clear
% close all
addpath(genpath('../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE

%% Load pre-computed FSP solutions + loaded data + MLEs:
% load('example_9_LoadingandFittingData_MLE.mat')

% View summary of 4-state STL1 model:
STL1_4state_MLE.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Metropolis-Hastings and Bayesian priors, iterating between 
%% computing the MLE and MH to sample uncertainty and improve model 
%% parameter fit to data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our 4-state STL1 model:
STL1_4state_MH = STL1_4state_MLE;

% Specify how many model parameters will be fit (the rest will be fixed):
fitpars = 13;

%% Specify Bayesian Prior and fit
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [0.5,2,5,3.5,-0.4,1,0.2,0.4,-0.5,-1.3,-0.1,2,0.5]; 

% Prior log-standard deviation:
sig_log10 = 2*ones(1,fitpars);  

% Prior:
STL1_4state_MH.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose parameters to search:
STL1_4state_MH.fittingOptions.modelVarsToFit = [1:fitpars]; 

% Create first parameter guess:
STL1_4state_MH_pars = [STL1_4state_MH.parameters{:,2}];      

% % You may need to re-run this multiple times until converged.
% % I got a MLE of -34,003.3 after a few runs. 

%% Iterating between MLE and MH
%  Running a few rounds of MLE and MH together may improve convergence.

STL1_4state_MH.parameters(:,2) = num2cell(STL1_4state_MH_pars);
for i=1:3
    % Maximize likelihood:
    STL1_4state_MH_pars = STL1_4state_MH.maximizeLikelihood([]);    
    % Update parameters in the model:
    STL1_4state_MH.parameters(1:fitpars, 2) = num2cell(STL1_4state_MH_pars);

    % Run Metropolis-Hastings    
    proposalWidthScale = 0.02;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

    % Set MH runtime options (number of samples, burnin, thin, etc.):
    MHOptions.numberOfSamples = 2000;
    MHOptions.burnin = 500;
    MHOptions.thin = 2;

    % Run Metropolis-Hastings (seeking acceptance ratio around 0.3-0.4): 
    [STL1_4state_MH_pars,~,STL1_4state_MHResults] = ...
        STL1_4state_MH.maximizeLikelihood([], MHOptions,...
        'MetropolisHastings');
    
    % Store MH parameters in model:
    STL1_4state_MH.parameters([1:fitpars],2) = num2cell(STL1_4state_MH_pars);
end
STL1_4state_MH.plotMHResults(STL1_4state_MHResults);

STL1_4state_MH.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1 (MH)', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=18, ProbXLim = [0 80],...
    TimePoints=[0 8 10 15 30 55], TitleFontSize=24, AxisLabelSize=20);

%% Save model & MH results:
saveNames = unique({ ...
    'STL1_4state_MH'
    'STL1_4state_MH_pars'
    'STL1_4state_MHResults'
    });
    
save('example_10_LoadingandFittingData_MH',saveNames{:})
