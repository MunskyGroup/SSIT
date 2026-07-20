%% SSIT/Examples/example_11b_LoadingandFittingData_MH_with_FIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.3.4: Loading and fitting time-varying STL1 yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%   * Use Bayesian priors and iterate between computing MLE and MH
%   * Use FIM for Metropolis-Hastings proposal distribution. 
%     This sometimes provides faster mixing (convergence).
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

% Make a new copy of our 4-state STL1 model:
STL1_4state_MH_FIM = STL1_4state;

%% Compute FIM, Run Metropolis Hastings
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [0.5,2,5,3.5,-0.4,1,0.2,0.4,-0.5,-1.3,-0.1,2,0.5]; 

% Prior log-standard deviation:
sig_log10 = 2*ones(1,13);      

% Prior:
STL1_4state_MH_FIM.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Compute total FIM including effect of prior:
fimTotal = STL1_4state_MH_FIM.evaluateExperiment(fimTotal_free,...
           STL1_4state_MH_FIM.dataSet.nCells,diag(sig_log10.^2)); 

% Select FIM for free parameters:
FIMfree = fimTotal{1}([1:13],[1:13]); 

% Estimate the covariance using CRLB:
COVfree = (1/2*(FIMfree + FIMfree'))^(-1);  

% Define Metropolis-Hasting settings:
STL1_4state_MH_FIM.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10([1:13])).^2./(2*sig_log10([1:13]).^2));
proposalWidthScale = 0.01;
STL1_4state_MH_FIM_FIMOptions = ...
 struct('proposalDistribution',@(x)mvnrnd(x,proposalWidthScale*COVfree),...
        'numberOfSamples',2000,'burnin',500,'thin',2);
fitOptions = optimset('Display','iter','MaxIter',2000);

% Compute the MLEs:
STL1_4state_MH_FIM = ...
    STL1_4state_MH_FIM.maximizeLikelihood(fitOptions=fitOptions);

% Run Metropolis-Hastings (seeking acceptance ratio around 0.3-0.4): 
STL1_4state_MH_FIM = STL1_4state_MH_FIM.maximizeLikelihood(...
    fitOptions=STL1_4state_MH_FIM_FIMOptions,...
    fitAlgorithm='MetropolisHastings');

% Plot MH samples, FIM:
STL1_4state_MH_FIM.plotMHResults(STL1_4state_MH_FIM.Solutions.mhResults,...
                                  FIM=FIMfree, fimScale='log')

STL1_4state_MH_FIM.plotFits(plotType="all",lineProps={'linewidth',2},...
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
