%% SSIT/Examples/example_10_LoadingandFittingData_MH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.3.4: Loading and fitting time-varying STL1 yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, data loaded in 
% example_8_LoadingandFittingData_DataLoading, and MLE computed in
% example_9_LoadingandFittingData_MLE
% clear
% close all

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE

%% Load pre-computed FSP solutions + loaded data + MLEs:
load('example_9_LoadingandFittingData_MLE.mat')

% View summary of 4-state STL1 model:
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Metropolis-Hastings and Bayesian priors, iterating between 
%% computing the MLE and MH to sample uncertainty and improve model 
%% parameter fit to data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify how many model parameters will be fit (the rest will be fixed):
fitpars = 13;
STL1_4state.fittingOptions.modelVarsToFit = [1:fitpars]; 

%% Specify Bayesian Prior and fit
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [0.5,2,5,3.5,-0.4,1,0.2,0.4,-0.5,-1.3,-0.1,2,0.5]; 

% Prior log-standard deviation:
sig_log10 = 2*ones(1,fitpars);  

% Prior:
STL1_4state.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));   

%% Iterating between MLE and MH
%  Running a few rounds of MLE and MH together may improve convergence.

% Set fitOptions, with the maximum allowable number of iterations to fit:
fitOptions = optimset('Display','iter','MaxIter',2000);

for i=1:2
    % Maximize likelihood:
    STL1_4state = STL1_4state.maximizeLikelihood(fitOptions=fitOptions);    

    % Run Metropolis-Hastings    
    proposalWidthScale = 0.01;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

    % Set MH runtime options (number of samples, burnin, thin, etc.):
    MHOptions.numberOfSamples = 2000;
    MHOptions.burnin = 500;
    MHOptions.thin = 2;

    % Run Metropolis-Hastings (seeking acceptance ratio around 0.3-0.4): 
    STL1_4state = STL1_4state.maximizeLikelihood(fitOptions=MHOptions,...
        fitAlgorithm='MetropolisHastings');
end

% Plot results:
STL1_4state.plotMHResults(STL1_4state.Solutions.mhResults);

STL1_4state.plotFits(plotType="all",lineProps={'linewidth',2},...
    Title='4-state STL1 (MH)', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=18, ProbXLim = [0 80],...
    TimePoints=[0 8 10 15 30 55], TitleFontSize=24, AxisLabelSize=20);

%% Save model & MH results:
saveNames = unique({ ...
    'STL1_4state'
    });
    
save('example_10_LoadingandFittingData_MH',saveNames{:})
