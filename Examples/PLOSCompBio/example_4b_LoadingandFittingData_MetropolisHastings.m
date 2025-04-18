%% example_7_MetropolisHastings
% Example script to demonstrate how to use Metropolis-Hastings to sample
% uncertainty (and improve model parameter fit to data)
addpath(genpath('../../src'));

%% Preliminaries
% Load our models from example_1_CreateSSITModels; compute FSP solutions 
% using example_2_SolveSSITModels_FSP and FIMs using example_4_FIM
%% Comment out the following 3 lines if example_4_FIM has already been run:
% clear
% close all
% example_3_FIM

% View model summaries
Model_FIM.summarizeModel
STL1_FIM.summarizeModel

% Make new copies of our models
Model_MH = Model_FIM;
STL1_MH = STL1_FIM;

%% Model: FIM inverse
% The inverse of the FIM provides an estimate of the model uncertainty.
% Here we are going to look at the FIM for the log of the model parameters
% and use that to compute the covariance of the log of the parameters.
% (Because parameters are positive values, but can very significantly in
% their magnitudes, it is often useful to examine them in a log-scale).
Model_FIMlog = diag([Model_MH.parameters{:,2}]) * Model_fimTotal{1}...
                                    * diag([Model_MH.parameters{:,2}]);
Model_covLog = Model_FIMlog^-1;

% The eigenvalues of covLog tells us what to expect for the uncertainty in
% the parameters.
[Model_eigVec,Model_eigVal] = eig(Model_covLog);
Model_eigVal = diag(Model_eigVal)

% Here, we see that there is one large direction of uncertainty, but the
% rest are pretty well constrained.  The direction of the greatest
% uncertainty is:
[~,j] = max(Model_eigVal);
Model_largestEigVec = Model_eigVec(:,j)

%% Model: Metropolis Hastings
% Now that we have an estimate of the shape of the uncertainty using the
% FIM we can now search parameter space and see what other parameter
% combinations are also closely matching to our data.  For this, we are
% going to use the Metropolis Hastings algorithm, where the proposal
% distribution is a multi-variate gaussian with a covariance that is
% proportional to the inverse FIM.  

% But because the FIM has some very small eigenvalues, we better may be
% better off reduceing the step size in those directions.  Here, we set it
% to at most one order of magnitude by adding an identity matrix to the FIM
% before inverting.
Model_covLogMod = (Model_FIMlog+1*diag(size(Model_FIMlog,1)))^(-1);

% Here, we set up the MH parameters:
Model_MH.solutionScheme = 'FSP'; % Set solutions scheme to FSP Sensitivity
Model_MH.fittingOptions.modelVarsToFit = 1:4;
Model_MHOptions = struct('numberOfSamples',3000,'burnin',100,'thin',3);
proposalWidthScale = 0.0000000000000001;
Model_MHOptions.proposalDistribution  = ...
 @(x)mvnrnd(x,proposalWidthScale * (Model_covLogMod + Model_covLogMod')/2);

% Next, we call the codes to sample the posterior parameter space:
[Modelpars,Model_likelihood,Model_chainResults] = ...
    Model_MH.maximizeLikelihood([],Model_MHOptions,'MetropolisHastings');
% When this runs, you want to see an acceptance of about 0.3 to 0.4, meaning that
% about a third of the proposals are accepted.  If the number is too small
% you need to decrease the proposal width; if it is too large you may need
% to increase the proposal width. For the default data set and model, I
% found that a scale of .5 to 5% of the FIM-based COV led to an okay acceptance
% rate, but this is variable and will change depending on the initial value
% in the chain.

% And now to plot the MH results and compare to the FIM.
Model_MH.plotMHResults(Model_chainResults,Model_FIMlog);

% Often the MH search can reveal a better parameter set, so let's make sure
% to update our model if it does:
Model_MH.parameters(:,2) = num2cell(Modelpars);
Model_MH.makeFitPlot
% If you do notice better fits, it would be good to re-run the fminsearch
% again - it is possible you can still find a better model to explain your
% data.  This can take several rounds of iteration before convergence.  I
% recommend creating a while loop to make it automated.

%% Model: Iterating between MLE and MH
% Let's run a few rounds of MLE and MH to see if we can get better
% convergence.
Model_MH.parameters(:,2) = num2cell(Modelpars);
for i=1:3
    % Maximize likelihood
    Modelpars = Model_MH.maximizeLikelihood([],fitOptions);    
    % Update parameters in the model:
    Model_MH.parameters(:,2) = num2cell(Modelpars);
    
    % Compute FIM
    Model_MH.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
    [Model_sensSoln] = Model_MH.solve;  % Solve the sensitivity problem
    Model_fimResults = Model_MH.computeFIM(Model_sensSoln.sens,'log');
    Model_FIMlog = Model_FIM.evaluateExperiment(Model_fimResults,Model_MH.dataSet.nCells);

    % Run Met. Hast.    
    Model_covLogMod = (Model_FIMlog{1} + diag(size(Model_FIMlog{1},1)))^(-1); % Adjusted proposal dist. covariance.
    proposalWidthScale = 0.0000000000000001;
    Model_MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(Model_covLogMod+Model_covLogMod')/2);
    [Modelpars,Model_likelihood,Model_chainResults] = Model_MH.maximizeLikelihood([],Model_MHOptions,'MetropolisHastings');
    % Update parameters in the model:
    Model_MH.parameters(:,2) = num2cell(Modelpars);
end
Model_MH.plotMHResults(Model_chainResults,Model_FIMlog);
Model_MH.makeFitPlot

%% Model: Evaluating the MH results
% Here we will generate three plots.  The first one will show the
% Likelihood function as we search over parameter space.  In order to get a
% good estimate of the parameter uncertainty, we want this to quickly reach
% ts maximum value and then to fluctuate around that value for a
% significant amount of time.  If you see that it is still incereasing, you
% know that the fit has not yet converged.
figure
subplot(1,3,1)
plot(Model_chainResults.mhValue)
title('MH Convergence')
xlabel('Iteration Number')
ylabel('LogLikelihood')

%% Compute FIM
Model_MH.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[Model_sensSoln] = Model_MH.solve;  % Solve the sensitivity problem
Model_fimResultsLog = Model_MH.computeFIM(Model_sensSoln.sens,'log');
Model_FIMlog = Model_MH.evaluateExperiment(Model_fimResultsLog,Model_MH.dataSet.nCells);
Model_fimResults= Model_MH.computeFIM(Model_sensSoln.sens,'lin');
Model_FIM = Model_MH.evaluateExperiment(Model_fimResults,Model_MH.dataSet.nCells);

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.

% Choose which parameters to compare.
Q = [3,4];  
subplot(1,3,2)

% Plot uncertainty in linear scale
Model_MH.makeMleFimPlot(exp(Model_chainResults.mhSamples)',Model_FIM{1},Q,0.95,1); hold on
title('Posterior -- Linear Scale')
xlabel(Model_MH.parameters{Q(1),1})
ylabel(Model_MH.parameters{Q(2),1})

% Plot uncertainty in log scale
subplot(1,3,3)
Model_MH.makeMleFimPlot(Model_chainResults.mhSamples',Model_FIMlog{1},Q,0.95,1); hold on
title('Posterior -- Natural Log Scale')
xlabel(['log ',Model_MH.parameters{Q(1),1}])
ylabel(['log ',Model_MH.parameters{Q(2),1}]) 

%% Model: Effective Sample Size
% For the MH, it is important to get a sense of how well it has sampled the
% posterior.  For this, we determine the effecitve sample size (i.e., the
% nuBmber of effectively indpendent samples within the MH chain).  This is
% found by examining the autocorrelation of the parameter chain to figure
% out the number of steps needed for correlations to decay and then divide
% the total number of steps by the de-correlation step.
figure
ipar = 4;
ac = xcorr(Model_chainResults.mhSamples(:,ipar)-mean(Model_chainResults.mhSamples(:,ipar)),'normalized');
ac = ac(size(Model_chainResults.mhSamples,1):end);
plot(ac,'LineWidth',3)
N = size(Model_chainResults.mhSamples,1);
tau = 1+2*sum(abs(ac(2:N/5)));
Neff = N/tau


%% STL1: FIM inverse
% The inverse of the FIM provides an estimate of the model uncertainty.
% Here we are going to look at the FIM for the log of the model parameters
% and use that to compute the covariance of the log of the parameters.
% (Because parameters are positive values, but can very significantly in
% their magnitudes, it is often useful to examine them in a log-scale).
STL1_FIMlog = diag([STL1_MH.parameters{:,2}]) * STL1_fimTotal{1}...
                    * diag([STL1_MH.parameters{:,2}]);
STL1_covLog = STL1_FIMlog^-1;

% The eigenvalues of covLog tells us what to expect for the uncertainty in
% the parameters.
[STL1_eigVec,STL1_eigVal] = eig(STL1_covLog);
STL1_eigVal = diag(STL1_eigVal)

% Here, we see that there is one large direction of uncertainty, but the
% rest are pretty well constrained.  The direction of the greatest
% uncertainty is:
[~,j] = max(STL1_eigVal);
STL1_largestEigVec = STL1_eigVec(:,j)

%% STL1: Metropolis Hastings
% Now that we have an estimate of the shape of the uncertainty using the
% FIM we can now search parameter space and see what other parameter
% combinations are also closely matching to our data.  For this, we are
% going to use the Metropolis Hastings algorithm, where the proposal
% distribution is a multi-variate gaussian with a covariance that is
% proportional to the inverse FIM.  

% But because the FIM has some very small eigenvalues, we better may be
% better off reduceing the step size in those directions.  Here, we set it
% to at most one order of magnitude by adding an identity matrix to the FIM
% before inverting.
STL1_covLogMod = (STL1_FIMlog + 1 * diag(size(STL1_FIMlog,1)))^(-1);

% Here, we set up the MH parameters:
STL1_MH.solutionScheme = 'FSP'; % Set solutions scheme to FSP Sensitivity
STL1_MH.fittingOptions.modelVarsToFit = 1:7;
STL1_MHOptions = struct('numberOfSamples',1000,'burnin',10,'thin',3);
proposalWidthScale = 0.00000001;
STL1_MHOptions.proposalDistribution = ...
  @(x)mvnrnd(x,proposalWidthScale * (STL1_covLogMod + STL1_covLogMod')/2);

% Next, we call the codes to sample the posterior parameter space:
[STL1pars,STL1_likelihood,STL1_chainResults] = ...
    STL1_MH.maximizeLikelihood([],STL1_MHOptions,'MetropolisHastings');
% When this runs, you want to see an acceptance of about 0.3 to 0.4, meaning that
% about a third of the proposals are accepted.  If the number is too small
% you need to decrease the proposal width; if it is too large you may need
% to increase the proposal width. For the default data set and model, I
% found that a scale of .5 to 5% of the FIM-based COV led to an okay acceptance
% rate, but this is variable and will change depending on the initial value
% in the chain.

% And now to plot the MH results and compare to the FIM.
STL1_MH.plotMHResults(STL1_chainResults,STL1_FIMlog);

% Often the MH search can reveal a better parameter set, so let's make sure
% to update our model if it does:
STL1_MH.parameters(:,2) = num2cell(STL1pars);
STL1_MH.makeFitPlot
% If you do notice better fits, it would be good to re-run the fminsearch
% again - it is possible you can still find a better model to explain your
% data.  This can take several rounds of iteration before convergence.  I
% recommend creating a while loop to make it automated.

%% STL1: Iterating between MLE and MH
% Let's run a few rounds of MLE and MH to see if we can get better
% convergence.
STL1_MH.parameters(:,2) = num2cell(STL1pars);
for i=1:3
    % Maximize likelihood
    STL1pars = STL1_MH.maximizeLikelihood([],fitOptions);    
    % Update parameters in the model:
    STL1_MH.parameters(:,2) = num2cell(STL1pars);
    
    % Compute FIM
    STL1_MH.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
    [STL1_sensSoln] = STL1_MH.solve;  % Solve the sensitivity problem
    STL1_fimResults = STL1_MH.computeFIM(STL1_sensSoln.sens,'log');
    STL1_FIMlog = STL1_MH.evaluateExperiment(STL1_fimResults,STL1_MH.dataSet.nCells);

    % Run Met. Hast.    
    STL1_covLogMod = (STL1_FIMlog{1} + diag(size(STL1_FIMlog{1},1)))^(-1); % Adjusted proposal dist. covariance.
    proposalWidthScale = 0.00000001;
    STL1_MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(STL1_covLogMod+STL1_covLogMod')/2);
    [STL1pars,STL1_likelihood,STL1_chainResults] = STL1_MH.maximizeLikelihood([],STL1_MHOptions,'MetropolisHastings');
    % Update parameters in the model:
    STL1_MH.parameters(:,2) = num2cell(STL1pars);
end
STL1_MH.plotMHResults(STL1_chainResults,STL1_FIMlog);
STL1_MH.makeFitPlot

%% STL1: Evaluating the MH results
% Here we will generate three plots.  The first one will show the
% Likelihood function as we search over parameter space.  In order to get a
% good estimate of the parameter uncertainty, we want this to quickly reach
% ts maximum value and then to fluctuate around that value for a
% significant amount of time.  If you see that it is still incereasing, you
% know that the fit has not yet converged.
figure
subplot(1,3,1)
plot(STL1_chainResults.mhValue)
title('MH Convergence')
xlabel('Iteration Number')
ylabel('LogLikelihood')

% Compute FIM
STL1_MH.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[STL1_sensSoln] = STL1_MH.solve;  % Solve the sensitivity problem
STL1_fimResultsLog = STL1_MH.computeFIM(STL1_sensSoln.sens,'log');
STL1_FIMlog = STL1_MH.evaluateExperiment(STL1_fimResultsLog,STL1_MH.dataSet.nCells);
STL1_fimResults= STL1_MH.computeFIM(STL1_sensSoln.sens,'lin');
STL1_FIM = STL1_MH.evaluateExperiment(STL1_fimResults,STL1_MH.dataSet.nCells);

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.

% Choose which parameters to compare.
Q = [3,4];  
subplot(1,3,2)

% Plot uncertainty in linear scale
STL1_MH.makeMleFimPlot(exp(STL1_chainResults.mhSamples)',STL1_FIM{1},Q,0.95,1); hold on
title('Posterior -- Linear Scale')
xlabel(STL1_MH.parameters{Q(1),1})
ylabel(STL1_MH.parameters{Q(2),1})

% Plot uncertainty in log scale
subplot(1,3,3)
STL1_MH.makeMleFimPlot(STL1_chainResults.mhSamples',STL1_FIMlog{1},Q,0.95,1); hold on
title('Posterior -- Natural Log Scale')
xlabel(['log ',STL1_MH.parameters{Q(1),1}])
ylabel(['log ',STL1_MH.parameters{Q(2),1}]) 

%% STL1: Effective Sample Size
% For the MH, it is important to get a sense of how well it has sampled the
% posterior.  For this, we determine the effecitve sample size (i.e., the
% nuBmber of effectively indpendent samples within the MH chain).  This is
% found by examining the autocorrelation of the parameter chain to figure
% out the number of steps needed for correlations to decay and then divide
% the total number of steps by the de-correlation step.
figure
ipar = 5;
ac = xcorr(STL1_chainResults.mhSamples(:,ipar)-mean(STL1_chainResults.mhSamples(:,ipar)),'normalized');
ac = ac(size(STL1_chainResults.mhSamples,1):end);
plot(ac,'LineWidth',3)
N = size(STL1_chainResults.mhSamples,1);
tau = 1+2*sum(abs(ac(2:N/5)));
Neff = N/tau