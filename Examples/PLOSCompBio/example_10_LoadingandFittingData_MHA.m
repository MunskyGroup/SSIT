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
example_9_LoadingandFittingData_MLE

% View model summaries:
STL1_MLE_refit.summarizeModel
STL1_MLE_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Metropolis-Hastings to sample uncertainty 
%   (and improve model parameter fit to data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of our STL1 model for Metropolis-Hastings (MH):
STL1_MH = STL1_MLE_refit;
STL1_4state_MH = STL1_MLE_4state;

    % Adjust proposal width scale (the default proposal distribution in 
    % SSIT is "@(x)x+0.1*randn(size(x))", which leads to low acceptance in
    % this case. A rule of thumb, depending on the data being analyzed, is 
    % to aim for an MH acceptance ratio of around 0.3-0.4 (some say 25%):
    proposalWidthScale = 0.01;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

% Set MH runtime options (number of samples, burnin, thin, etc.):
MHOptions.numberOfSamples = 3000;
MHOptions.burnin = 300;
MHOptions.thin = 3;

% Run Metropolis-Hastings: 
[STL1_MH_pars,~,STL1_MHResults] = STL1_MH.maximizeLikelihood(...
    [], MHOptions, 'MetropolisHastings');

[STL1_MH_pars_4state,~,STL1_MHResults_4state] = ...
    STL1_4state_MH.maximizeLikelihood([], MHOptions, 'MetropolisHastings');

% Store MH parameters in model:
STL1_MH.parameters([1:7],2) = num2cell(STL1_MH_pars);
STL1_4state_MH.parameters([1:15],2) = num2cell(STL1_MH_pars_4state);

% Plot results:
STL1_MH.plotMHResults(STL1_MHResults,[],'log',[])
STL1_MH.makeFitPlot

STL1_4state_MH.plotMHResults(STL1_MHResults_4state,[],'log',[])
STL1_4state_MH.makeFitPlot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Use FIM for Metropolis-Hastings proposal distribution 
%   This sometimes provides faster mixing (convergence), although in our 
%   simple example, the default proposal distribution (above) is fine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make new copies of our models:
STL1_MH_FIM = STL1_MLE_refit;

%% Compute FIM, Run Metropolis Hastings
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [-1,1,1,1,1,-2,-2];  

% Prior log-standard deviation:
sig_log10 = 2*ones(1,7);      

% Prior:
STL1_MH_FIM.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose parameters to search:
STL1_MH_FIM.fittingOptions.modelVarsToFit = [1:7];

% Create first parameter guess:
STL1_MH_FIM_pars = [STL1_MH_FIM.parameters{:,2}];         

% Compute individual FIMs:
fimResults = STL1_MH_FIM.computeFIM([],'log'); 

% Compute total FIM including effect of prior:
fimTotal = STL1_MH_FIM.evaluateExperiment(fimResults,...
    STL1_MH_FIM.dataSet.nCells,diag(sig_log10.^2)); 

% Choose parameters to search:
STL1_MH_FIM.fittingOptions.modelVarsToFit = [1:7]; 

% Select FIM for free parameters:
FIMfree = fimTotal{1}([1:7],[1:7]); 

% Estimate the covariance using CRLB:
COVfree = (1/2*(FIMfree+FIMfree'))^(-1);  

% Define Metropolis-Hasting settings:
STL1_MH_FIM.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10([1:7])).^2./(2*sig_log10([1:7]).^2));
proposalWidthScale = 0.0001;
STL1_MH_FIMOptions = ...
 struct('proposalDistribution',@(x)mvnrnd(x,proposalWidthScale*COVfree),...
 'numberOfSamples',3000,'burnin',300,'thin',3);

% Run Metropolis Hastings
[STL1_MH_FIM_pars,~,MHResultsDusp1] = STL1_MH_FIM.maximizeLikelihood(...
        [], STL1_MH_FIMOptions, 'MetropolisHastings'); 

% Store sampled parameters:
STL1_MH_FIM.parameters([1:7],2) = num2cell(STL1_MH_FIM_pars);

% Plot MH samples, FIM:
STL1_MH_FIM.plotMHResults(MHResultsDusp1,FIMfree,'log',[])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Use Bayesian priors and iterate between computing MLE and MH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify Bayesian Prior and fit
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [-1,1,1,1,1,-2,-2]; 

% Prior log-standard deviation:
sig_log10 = 2*ones(1,7);  

% Prior:
STL1_MH_FIM.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose parameters to search:
STL1_MH_FIM.fittingOptions.modelVarsToFit = [1:7]; 

% Create first parameter guess:
STL1_MH_FIM_pars = [STL1_MH_FIM.parameters{:,2}];      

% Fit to maximize likelihood:
STL1_MH_FIM_pars = STL1_MH_FIM.maximizeLikelihood(STL1_MH_FIM_pars); 

% Update new parameters:
STL1_MH_FIM.parameters(:,2) = num2cell(STL1_MH_FIM_pars); 

% Plot fitting results:
STL1_MH_FIM.makeFitPlot  
% You may need to re-run this multiple times until converged.
% I got a MLE of -52,454.1 after a few runs. 


%% Iterating between MLE and MH
%  Running a few rounds of MLE and MH together may improve convergence.

STL1_MH_FIM.parameters(:,2) = num2cell(STL1_MH_FIM_pars);
for i=1:3
    % Maximize likelihood:
    STL1_MH_FIM_pars = STL1_MH_FIM.maximizeLikelihood([],fitOptions);    
    % Update parameters in the model:
    STL1_MH_FIM.parameters(:,2) = num2cell(STL1_MH_FIM_pars);

    % Compute FIM
    % Set solutions scheme to FSP Sensitivity:
    STL1_MH_FIM.solutionScheme = 'fspSens'; 
    % Solve the sensitivity problem:
    [STL1_sensSoln] = STL1_MH_FIM.solve;  
    STL1_fimResults = STL1_MH_FIM.computeFIM(STL1_sensSoln.sens,'log');
    STL1_FIMlog = STL1_MH_FIM.evaluateExperiment(STL1_fimResults,...
                                               STL1_MH_FIM.dataSet.nCells);

    % Run Metropolis-Hastings    
    % Adjusted proposal distribution covariance:
    STL1_covLogMod = (STL1_FIMlog{1} + diag(size(STL1_FIMlog{1},1)))^(-1); 
    proposalWidthScale = 0.00001;
    STL1_MH_FIMOptions.proposalDistribution  = ...
       @(x)mvnrnd(x,proposalWidthScale*(STL1_covLogMod+STL1_covLogMod')/2);
    [STL1_MH_FIM_pars,STL1_likelihood_FIM,STL1_chainResults_FIM] = ...
        STL1_MH_FIM.maximizeLikelihood([],STL1_MH_FIMOptions,...
                                        'MetropolisHastings');
    % Update parameters in the model:
    STL1_MH_FIM.parameters(:,2) = num2cell(STL1_MH_FIM_pars);
end
STL1_MH_FIM.plotMHResults(STL1_chainResults_FIM,STL1_FIMlog);
STL1_MH_FIM.makeFitPlot


%% Evaluating the MH results
% Here we will generate three plots.  The first one will show the
% likelihood function as we search over parameter space.  In order 
% to get a good estimate of the parameter uncertainty, we want this 
% to quickly reach ts maximum value and then to fluctuate around 
% that value for a significant amount of time.  If you see that it 
% is still incereasing, you know that the fit has not yet converged.
figure
plot(STL1_chainResults_FIM.mhValue)
title('MH Convergence')
xlabel('Iteration Number')
ylabel('LogLikelihood')

%% Compute FIM
% Set solutions scheme to FSP Sensitivity:
STL1_MH_FIM.solutionScheme = 'fspSens';
% Solve the sensitivity problem:
[STL1_sensSoln] = STL1_MH_FIM.solve;  
STL1_fimResultsLog = STL1_MH_FIM.computeFIM(STL1_sensSoln.sens,'log');
STL1_FIMlog = STL1_MH_FIM.evaluateExperiment(STL1_fimResultsLog,...
                                             STL1_MH_FIM.dataSet.nCells);
STL1_fimResults= STL1_MH_FIM.computeFIM(STL1_sensSoln.sens,'lin');
STL1_FIM = STL1_MH_FIM.evaluateExperiment(STL1_fimResults,...
                                          STL1_MH_FIM.dataSet.nCells);

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.

% Choose which parameters to compare:
Q = [3,4];  
subplot(1,3,2)

% Plot uncertainty in linear scale:
STL1_MH_FIM.makeMleFimPlot(exp(STL1_chainResults_FIM.mhSamples)',...
                             STL1_FIM{1},Q,0.95,1); hold on
title('Posterior -- Linear Scale')
xlabel(STL1_MH_FIM.parameters{Q(1),1})
ylabel(STL1_MH_FIM.parameters{Q(2),1})

% Plot uncertainty in log scale:
subplot(1,3,3)
STL1_MH_FIM.makeMleFimPlot(STL1_chainResults_FIM.mhSamples',...
                           STL1_FIMlog{1},Q,0.95,1); hold on
title('Posterior -- Natural Log Scale')
xlabel(['log ',STL1_MH_FIM.parameters{Q(1),1}])
ylabel(['log ',STL1_MH_FIM.parameters{Q(2),1}]) 

%% Effective Sample Size
% For the MH, it is important to get a sense of how well it has sampled the
% posterior.  For this, we determine the effective sample size (i.e., the
% number of effectively indpendent samples within the MH chain).  This is
% found by examining the autocorrelation of the parameter chain to figure
% out the number of steps needed for correlations to decay and then divide
% the total number of steps by the de-correlation step.
figure
ipar = 4;
ac = xcorr(STL1_chainResults_FIM.mhSamples(:,ipar) ...
     - mean(STL1_chainResults_FIM.mhSamples(:,ipar)),'normalized');
ac = ac(size(STL1_chainResults_FIM.mhSamples,1):end);
plot(ac,'LineWidth',3)
N = size(STL1_chainResults_FIM.mhSamples,1);
tau = 1+2*sum(abs(ac(2:N/5)));
Neff = N/tau