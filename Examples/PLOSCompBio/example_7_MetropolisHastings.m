%% example_7_MetropolisHastings
% Example script to demonstrate how to use Metropolis-Hastings to sample
% uncertainty
close all
addpath(genpath('../src'));

%% Preliminaries
% Load our models described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
example_6_LoadingandFittingData_MLE
Model.summarizeModel
STL1Model.summarizeModel

%%

[FIM] = Model_FIM.evaluateExperiment(fimResults,Model_FIM.dataSet.nCells);

% The inverse of the FIM provides an estimate of the model uncertainty.
% Here we are going to look at the FIM for the log of the model parameters
% and use that to compute the covariance of the log of the parameters.
% (Because parameters are positive values, but can very significantly in
% their magnitudes, it is often useful to examine them in a log-scale).
FIMlog = diag([Model_FIM.parameters{:,2}])*FIM{1}*diag([Model_FIM.parameters{:,2}]);
covLog = FIMlog^-1;

% The eigenvalues of covLog tells us what to expect for the uncertainty in
% the parameters.
[eigVec,eigVal] = eig(covLog);
eigVal = diag(eigVal)

% Here, we see that there is one large direction of uncertainty, but the
% rest are pretty well constrained.  The direction of the greatest
% uncertainty is:
[~,j] = max(eigVal);
largestEigVec = eigVec(:,j)

%% Metropolis Hastings
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
covLogMod = (FIMlog+1*diag(size(FIMlog,1)))^(-1);

% Here, we set up the MH parameters:
Model_FIM.solutionScheme = 'FSP'; % Set solutions scheme to FSP Sensitivity
Model_FIM.fittingOptions.modelVarsToFit = [1:4];
MHOptions = struct('numberOfSamples',1000,'burnin',0,'thin',3);
proposalWidthScale = 0.4;
MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(covLogMod+covLogMod')/2);

% Next, we call the codes to sample the posterior parameter space:
[Modelpars,likelihood,chainResults] = Model_FIM.maximizeLikelihood([],MHOptions,'MetropolisHastings');
% When this runs, you want to see an acceptance of about 0.3 to 0.4, meaning that
% about a third of the proposals are accepted.  If the number is too small
% you need to decrease the proposal width; if it is too large you may need
% to increase the proposal width. For the default data set and model, I
% found that a scale of .5 to 5% of the FIM-based COV led to an okay acceptance
% rate, but this is variable and will change depending on the initial value
% in the chain.

% And now to plot the MH results and compare to the FIM.
Model_FIM.plotMHResults(chainResults,FIMlog);

% Often the MH search can reveal a better parameter set, so let's make sure
% to update our model if it does:
Model_FIM.parameters(:,2) = num2cell(Modelpars);
Model_FIM.makeFitPlot
% If you do notice better fits, it would be good to re-run the fminsearch
% again - it is possible you can still find a better model to explain your
% data.  This can take several rounds of iteration before convergence.  I
% recommend creating a while loop to make it automated.

%% Iterating between MLE and MH.
% Let's run a few rounds of MLE and MH to see if we can get better
% convergence.
Model_FIM.parameters(:,2) = num2cell(Modelpars);
for i=1:3
    % Maximize likelihood
    pars = Model_FIM.maximizeLikelihood([],fitOptions);    
    % Update parameters in the model:
    Model_FIM.parameters(:,2) = num2cell(pars);
    
    % Compute FIM
    Model_FIM.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
    [sensSoln] = Model_FIM.solve;  % Solve the sensitivity problem
    fimResults = Model_FIM.computeFIM(sensSoln.sens,'log');
    FIMlog = Model_FIM.evaluateExperiment(fimResults,Model_FIM.dataSet.nCells);

    % Run Met. Hast.    
    covLogMod = (FIMlog{1} + diag(size(FIMlog{1},1)))^(-1); % Adjusted proposal dist. covariance.
    proposalWidthScale = 0.4;
    MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(covLogMod+covLogMod')/2);
    [pars,likelihood,chainResults] = Model_FIM.maximizeLikelihood([],MHOptions,'MetropolisHastings');
    % Update parameters in the model:
    Model_FIM.parameters(:,2) = num2cell(pars);
end
Model_FIM.plotMHResults(chainResults,FIMlog);
Model_FIM.makeFitPlot

%% Evaluating the MH results
% Here we will generate three plots.  The first one will show the
% Likelihood function as we search over parameter space.  In order to get a
% good estimate of the parameter uncertainty, we want this to quickly reach
% ts maximum value and then to fluctuate around that value for a
% significant amount of time.  If you see that it is still incereasing, you
% know that the fit has not yet converged.
figure
subplot(1,3,1)
plot(chainResults.mhValue)
title('MH Convergence')
xlabel('Iteration Number')
ylabel('LogLikelihood')

% Compute FIM
Model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln] = Model.solve;  % Solve the sensitivity problem
fimResultsLog = Model.computeFIM(sensSoln.sens,'log');
FIMlog = Model.evaluateExperiment(fimResultsLog,Model.dataSet.nCells);
fimResults= Model.computeFIM(sensSoln.sens,'lin');
FIM = Model.evaluateExperiment(fimResults,Model.dataSet.nCells);

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.

% Choose which parameters to compare.
Q = [3,4];  
subplot(1,3,2)

% Plot uncertainty in linear scale
Model.makeMleFimPlot(exp(chainResults.mhSamples)',FIM{1},Q,0.95,1); hold on
title('Posterior -- Linear Scale')
xlabel(Model.parameters{Q(1),1})
ylabel(Model.parameters{Q(2),1})

% Plot uncertainty in log scale
subplot(1,3,3)
Model.makeMleFimPlot(chainResults.mhSamples',FIMlog{1},Q,0.95,1); hold on
title('Posterior -- Natural Log Scale')
xlabel(['log ',Model.parameters{Q(1),1}])
ylabel(['log ',Model.parameters{Q(2),1}]) 

%% Effective Sample Size
% For the MH, it is important to get a sense of how well it has sampled the
% posterior.  For this, we determine the effecitve sample size (i.e., the
% nuBmber of effectively indpendent samples within the MH chain).  This is
% found by examining the autocorrelation of the parameter chain to figure
% out the number of steps needed for correlations to decay and then divide
% the total number of steps by the de-correlation step.
figure
ipar = 5;
ac = xcorr(chainResults.mhSamples(:,ipar)-mean(chainResults.mhSamples(:,ipar)),'normalized');
ac = ac(size(chainResults.mhSamples,1):end);
plot(ac,'LineWidth',3)
N = size(chainResults.mhSamples,1);
tau = 1+2*sum(abs(ac(2:N/5)));
Neff = N/tau


%%

[FIM] = STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_FIM.dataSet.nCells);

% The inverse of the FIM provides an estimate of the model uncertainty.
% Here we are going to look at the FIM for the log of the model parameters
% and use that to compute the covariance of the log of the parameters.
% (Because parameters are positive values, but can very significantly in
% their magnitudes, it is often useful to examine them in a log-scale).
FIMlog = diag([STL1_FIM.parameters{:,2}])*FIM{1}*diag([STL1_FIM.parameters{:,2}]);
covLog = FIMlog^-1;

% The eigenvalues of covLog tells us what to expect for the uncertainty in
% the parameters.
[eigVec,eigVal] = eig(covLog);
eigVal = diag(eigVal)

% Here, we see that there is one large direction of uncertainty, but the
% rest are pretty well constrained.  The direction of the greatest
% uncertainty is:
[~,j] = max(eigVal);
largestEigVec = eigVec(:,j)

%% Metropolis Hastings
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
covLogMod = (FIMlog+1*diag(size(FIMlog,1)))^(-1);

% Here, we set up the MH parameters:
STL1_FIM.solutionScheme = 'FSP'; % Set solutions scheme to FSP Sensitivity
STL1_FIM.fittingOptions.modelVarsToFit = [1:8];
MHOptions = struct('numberOfSamples',1000,'burnin',0,'thin',3);
proposalWidthScale = 0.4;
MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(covLogMod+covLogMod')/2);

% Next, we call the codes to sample the posterior parameter space:
[STL1pars,STL1likelihood,STL1chainResults] = STL1_FIM.maximizeLikelihood([],MHOptions,'MetropolisHastings');
% When this runs, you want to see an acceptance of about 0.3 to 0.4, meaning that
% about a third of the proposals are accepted.  If the number is too small
% you need to decrease the proposal width; if it is too large you may need
% to increase the proposal width. For the default data set and model, I
% found that a scale of .5 to 5% of the FIM-based COV led to an okay acceptance
% rate, but this is variable and will change depending on the initial value
% in the chain.

% And now to plot the MH results and compare to the FIM.
STL1_FIM.plotMHResults(STL1chainResults,FIMlog);

% Often the MH search can reveal a better parameter set, so let's make sure
% to update our model if it does:
STL1_FIM.parameters(:,2) = num2cell(STL1pars);
STL1_FIM.makeFitPlot
% If you do notice better fits, it would be good to re-run the fminsearch
% again - it is possible you can still find a better model to explain your
% data.  This can take several rounds of iteration before convergence.  I
% recommend creating a while loop to make it automated.

%% Iterating between MLE and MH.
% Let's run a few rounds of MLE and MH to see if we can get better
% convergence.
STL1_FIM.parameters(:,2) = num2cell(STL1pars);
for i=1:3
    % Maximize likelihood
    STL1pars = STL1_FIM.maximizeLikelihood([],fitOptions);    
    % Update parameters in the model:
    STL1_FIM.parameters(:,2) = num2cell(STL1pars);
    
    % Compute FIM
    STL1_FIM.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
    [sensSoln] = STL1_FIM.solve;  % Solve the sensitivity problem
    STL1fimResults = STL1_FIM.computeFIM(sensSoln.sens,'log');
    FIMlog = STL1_FIM.evaluateExperiment(STL1fimResults,STL1_FIM.dataSet.nCells);

    % Run Met. Hast.    
    covLogMod = (FIMlog{1} + diag(size(FIMlog{1},1)))^(-1); % Adjusted proposal dist. covariance.
    proposalWidthScale = 0.4;
    MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(covLogMod+covLogMod')/2);
    [pars,likelihood,chainResults] = Model.maximizeLikelihood([],MHOptions,'MetropolisHastings');
    % Update parameters in the model:
    Model.parameters(:,2) = num2cell(pars);
end
Model.plotMHResults(chainResults,FIMlog);
Model.makeFitPlot

%% Evaluating the MH results
% Here we will generate three plots.  The first one will show the
% Likelihood function as we search over parameter space.  In order to get a
% good estimate of the parameter uncertainty, we want this to quickly reach
% ts maximum value and then to fluctuate around that value for a
% significant amount of time.  If you see that it is still incereasing, you
% know that the fit has not yet converged.
figure
subplot(1,3,1)
plot(chainResults.mhValue)
title('MH Convergence')
xlabel('Iteration Number')
ylabel('LogLikelihood')

% Compute FIM
Model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln] = Model.solve;  % Solve the sensitivity problem
fimResultsLog = Model.computeFIM(sensSoln.sens,'log');
FIMlog = Model.evaluateExperiment(fimResultsLog,Model.dataSet.nCells);
fimResults= Model.computeFIM(sensSoln.sens,'lin');
FIM = Model.evaluateExperiment(fimResults,Model.dataSet.nCells);

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.

% Choose which parameters to compare.
Q = [3,4];  
subplot(1,3,2)

% Plot uncertainty in linear scale
STL1_FIM.makeMleFimPlot(exp(chainResults.mhSamples)',FIM{1},Q,0.95,1); hold on
title('Posterior -- Linear Scale')
xlabel(Model.parameters{Q(1),1})
ylabel(Model.parameters{Q(2),1})

% Plot uncertainty in log scale
subplot(1,3,3)
Model.makeMleFimPlot(chainResults.mhSamples',FIMlog{1},Q,0.95,1); hold on
title('Posterior -- Natural Log Scale')
xlabel(['log ',Model.parameters{Q(1),1}])
ylabel(['log ',Model.parameters{Q(2),1}]) 

%% Effective Sample Size
% For the MH, it is important to get a sense of how well it has sampled the
% posterior.  For this, we determine the effecitve sample size (i.e., the
% nuBmber of effectively indpendent samples within the MH chain).  This is
% found by examining the autocorrelation of the parameter chain to figure
% out the number of steps needed for correlations to decay and then divide
% the total number of steps by the de-correlation step.
figure
ipar = 5;
ac = xcorr(chainResults.mhSamples(:,ipar)-mean(chainResults.mhSamples(:,ipar)),'normalized');
ac = ac(size(chainResults.mhSamples,1):end);
plot(ac,'LineWidth',3)
N = size(chainResults.mhSamples,1);
tau = 1+2*sum(abs(ac(2:N/5)));
Neff = N/tau