%% SSIT/Examples/example_17_ABC.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5.7: Approximate Bayesian Computation 
%
% Demonstration of Approximate Bayesian Computation (ABC) in the SSIT
% using the runABCsearch method.
%
% This example:
%   1. Loads our 4-state STL1 model with SSA solution scheme.
%   2. Associates the 4-state STL1 model with STL1 smFISH data.
%   3. Defines a prior over parameters.
%   4. Runs ABC via Metropolis–Hastings using 'cdf_one_norm' loss.
%   5. Visualizes the ABC results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% clear
% close all
addpath(genpath('../src'));

%% Set SSA options:
% Make copy of our 4-state STL1 model:
STL1_4state_ABC = STL1_4state_MH;

%Set solution scheme to SSA:
STL1_4state_ABC.solutionScheme = 'SSA';

% Set number of simulations performed per experiment (small # for demo):
STL1_4state_ABC.ssaOptions.nSimsPerExpt=100;
    
% Equilibrate before starting (burn-in):
STL1_4state_ABC.tSpan = [0,STL1_4state_ABC.tSpan];

% Run iterations in parallel with multiple cores:
STL1_4state_ABC.ssaOptions.useParallel = true;

%% Associate STL1 data:
STL1_4state_ABC = ...
    STL1_4state_ABC.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                            {'mRNA','RNA_STL1_total_TS3Full'},...
                            {'Replica',1;'Condition','0.2M_NaCl_Step'});

% Choose which parameters to fit (here: all parameters):
STL1_4state_ABC.fittingOptions.modelVarsToFit = ...
    1:size(STL1_4state_ABC.parameters,1);

%% Set up a prior over parameters (logPriorLoss)
% logPriorLoss should return a *loss* (positive penalty); smaller is better.
% A convenient choice is a quadratic penalty in log10-parameter space,
% corresponding to a log-normal prior.

% Get current parameter values as a reasonable prior mean:
theta0 = cell2mat(STL1_4state_ABC.parameters(:,2));
log10_mu = log10(theta0(:));
log10_sigma = 2 * ones(size(log10_mu));  % std dev in log10-space  

% Define prior "loss" (default, @(x)allFitOptions.obj(exp(x))):
logPriorLoss = [];

%% Set ABC / MCMC options
% runABCsearch passes 'fitOptions' to maximizeLikelihood with the
% 'MetropolisHastings' algorithm. Tune these depending on your problem size.

fitOptions = struct();
fitOptions.numberOfSamples       = 1000;         % Total MH iterations 
fitOptions.burnIn                = 100;          % Discard burn-in samples
fitOptions.thin                  = 0;           % Keep every nth sample
proposalWidthScale               = 0.2;        % Proposal scale
% Proposal distribution:
fitOptions.proposalDistribution  = @(x)x+proposalWidthScale*randn(size(x));
% Log prior:
fitOptions.logPrior = @(x)sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Initial parameter guess (optional, default: current Model.parameters):
parGuess = [];

% Choose loss function for ABC (default: 'cdf_one_norm'):
lossFunction = 'cdf_one_norm';

% Enforce independence by downsampling SSA trajectories:
enforceIndependence = true;

%% Run ABC search
% This will:
%   * repeatedly simulate SSA trajectories,
%   * compute a CDF-based loss against the data,
%   * add the prior penalty, and
%   * perform MH sampling to approximate the posterior.
%
% Outputs:
%   pars                    - “best” (minimum-loss) parameter set found
%   minimumLossFunction     - value of the loss at that point
%   Results                 - MH/ABC diagnostics and chains
%   ModelABC                - model updated with 'pars'

% Compile and store the given reaction propensities:
STL1_4state_ABC = STL1_4state_ABC.formPropensitiesGeneral('STL1_4state_ABC');

[parsABC, minimumLoss, ResultsABC, STL1_4state_ABC] = ...
     STL1_4state_ABC.runABCsearch(parGuess, lossFunction, logPriorLoss,...
                                  fitOptions, enforceIndependence);

fprintf('ABC completed.\n');
fprintf('Minimum loss value: %g\n', minimumLoss);
disp('Best-fit parameters (ABC):');
disp(parsABC(:).');

%% Inspect ABC results: 
% The 'ResultsABC' struct is returned by maximizeLikelihood with the
% 'MetropolisHastings' algorithm. 
%       ResultsABC.mhSamples        - MCMC chain of parameter samples
%       ResultsABC.mhValue          - corresponding loss values
%       ResultsABC.mhAcceptance     - MH acceptance fraction
% Below we show a simple marginal histogram for each fitted parameter.

if isfield(ResultsABC, 'mhSamples')
    parChain = ResultsABC.mhSamples;   % size: [numberOfSamples x nPars] 
    nPars    = size(parChain, 2);

    figure;
    for k = 1:nPars
        subplot(ceil(nPars/2), 2, k);
        histogram(parChain(:,k), 40, 'Normalization', 'pdf');
        hold on;
        xline(parsABC(k), 'r', 'LineWidth', 1.5);
        xline(STL1_4state_MH_pars(k), 'b', 'LineWidth', 1.5);
        title(sprintf('Parameter %d', k));
        xlabel('\theta_k');
        ylabel('Posterior density (approx.)');
    end
    sgtitle('ABC posterior marginals (approximate)');
else
    warning('ResultsABC.mhSamples not found.');
end

%% Compare initial vs. final (ABC) parameter losses 
%% (Experimental data vs. data simulated by the SSA):
% Minimum loss from ABC run:
minLoss = minimumLoss;  
nTimes  = sum(STL1_4state_ABC.fittingOptions.timesToFit);
% Set the number of species being fitting (in this case, only mRNA):
nSpecies = 1;  
% Replicate groups / dose groups etc., if applicable:
nConds   = 1; 

avgLossPerCDF = minLoss / (nTimes * nSpecies * nConds);
fprintf('Average CDF L1 discrepancy per time/species: %.4f\n',avgLossPerCDF);

L_init = STL1_4state_ABC.computeLossFunctionSSA(lossFunction,...
                                                theta0, enforceIndependence);
L_min  = minimumLoss;

fprintf('Initial loss: %.3f,  Final (min) loss: %.3f\n', L_init, L_min);
fprintf('Relative improvement: %.1f%%\n', 100 * (L_init - L_min)/L_init);


%% Compare ABC posterior sample to MLE
% TODO: Overlay parameter values and compute predictive distributions.

thetaMLE = cell2mat(STL1_4state_MLE.parameters(1:15,2));

L_MLE = STL1_4state_ABC.computeLossFunctionSSA(lossFunction,...
                                            thetaMLE, enforceIndependence);
L_min  = minimumLoss;

fprintf('MLE loss: %.3f,  Final (min) ABC loss: %.3f\n', L_MLE, L_min);
fprintf('Relative improvement: %.1f%%\n', 100 * (L_MLE - L_min)/L_MLE);


