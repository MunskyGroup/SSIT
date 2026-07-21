%% SSIT/Examples/example_11b_LoadingandFittingData_MH_with_FIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.3.4: Loading and fitting time-varying STL1 yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%   * Use Bayesian priors and iterate between computing MLE and MH
%   * Use FIM for Metropolis-Hastings proposal distribution. 
%     This often provides faster mixing (convergence).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, data loaded in 
% example_9_LoadingandFittingData_DataLoading, and MLE computed in
% example_10_LoadingandFittingData_MLE
% clear
% close all

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_7_FIM
% example_9_LoadingandFittingData_DataLoading
% example_10_LoadingandFittingData_MLE

%% Load FIM results and pre-computed FSP solutions + loaded data + MLEs:
load('example_7_FIM.mat')
load('example_10_LoadingandFittingData_MLE.mat')

%% Compute FIM, Run Metropolis Hastings
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [0.5,2,5,3.5,-0.4,1,0.2,0.4,-0.5,-1.3,-0.1,2,0.5]; 

% Prior log-standard deviation:
sig_log10 = 2*ones(1,13);      

% Prior:
STL1_4state.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Calculate the FIMs for all time points for the free parameters
fims_free = STL1_4state.computeFIM(scale='log',observed='mRNA',...
    freePars=freePars); % FIM sub matrix

% Compute total FIM including effect of prior. Note, this FIM calculation
% is  made at the current parameter values which were fit to maximize the
% MLE in example "example_10_LoadingandFittingData_MLE".
fimTotal = STL1_4state.evaluateExperiment(fims_free,...
           STL1_4state.dataSet.nCells,diag(sig_log10.^2)*log(10)^2); 

% Select FIM for free parameters:
FIMfree = fimTotal{1}; 

% Estimate the covariance using CRLB:
COVfree = (1/2*(FIMfree + FIMfree'))^(-1);  

% Define Metropolis-Hasting settings:
STL1_4state.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10((1:13))).^2./(2*sig_log10((1:13)).^2));
proposalWidthScale = 1;
STL1_4state_MH_FIM_FIMOptions = ...
 struct('proposalDistribution',@(x)mvnrnd(x,proposalWidthScale*COVfree),...
        'numberOfSamples',2000,'burnin',500,'thin',2);
fitOptions = optimset('Display','iter','MaxIter',2000);

% Run Metropolis-Hastings (seeking acceptance ratio around 0.3-0.4): 
STL1_4state = STL1_4state.maximizeLikelihood(...
    fitOptions=STL1_4state_MH_FIM_FIMOptions,...
    fitAlgorithm='MetropolisHastings');

% Plot MH samples and FIM:
STL1_4state.plotMHResults(STL1_4state.Solutions.mhResults,...
                                  FIM=FIMfree, fimScale='log')

STL1_4state.plotFits(plotType="all",lineProps={'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

%% Save models & MH results:
saveNames = unique({'STL1_4state_MH_FIM'
    'fimResults'
    'fimTotal'
    'FIMfree'
    'COVfree'
    });
    
save('example_11b_LoadingandFittingData_MH_with_FIM',saveNames{:})
