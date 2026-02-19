%% SSIT/Examples/example_10b_LoadingandFittingData_MH_with_FIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.3: Loading and fitting time-varying STL1 yeast data 
%   * Uncertainty sampling using the Metropolis-Hastings Algorithm (MHA)
%   * Use Bayesian priors and iterate between computing MLE and MH
%   * Use FIM for Metropolis-Hastings proposal distribution 
%     This sometimes provides faster mixing (convergence), although in our 
%     simple example, the default proposal distribution (above) is fine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our 4-state STL1 model:
STL1_4state_MH_FIM = STL1_4state_MLE;

%% Compute FIM, Run Metropolis Hastings
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [0.8,3,-0.1,2,2.75,0.6,3,2.5,0,3.5,1.5,-0.15,0.5,1.5,-1]; 

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
[STL1_4state_MH_FIM_pars,~,STL1_4state_FIM_MHResults] = ...
    STL1_4state_MH_FIM.maximizeLikelihood([], ...
    STL1_4state_MH_FIM_FIMOptions, 'MetropolisHastings'); 

% Store sampled parameters:
STL1_4state_MH_FIM.parameters([1:15],2) = ...
    num2cell(STL1_4state_MH_FIM_pars);

% Plot MH samples, FIM:
STL1_4state_MH_FIM.plotMHResults(STL1_4state_FIM_MHResults,FIMfree,'log',[])

STL1_4state_MH_FIM.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

%% Save models & MH results:
saveNames = unique({'STL1_4state_MH_FIM'
    'STL1_4state_MH_FIM_pars'
    'fimResults'
    'fimTotal'
    'FIMfree'
    'COVfree'
    'STL1_4state_FIM_MHResults'
    });
    
save('example_10b_LoadingandFittingData_MH_with_FIM',saveNames{:})
