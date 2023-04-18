%% Example script to show SSIT for the Vanderbilt Q_BIO Group
% In this script, we are going to show how to create, solve and fit a CME
% model to some single-cell smFISH data.  For this example, we will use
% some data collected in Dr. Gregor Neuert's laboratory at Vanderbilt.
clear all
close all
clc

%% Create SSIT Model
% First, we are going to create an FSP model for a bursting gene expression
% model.  This model will consist of 3 species: OFF Gene, ON Gene, and
% mRNA.  There will be four reactions: activation, inactivation,
% transcription and degradation.  The activation rate will be assumed to be
% time varying and controlled by a MAPK signal.

Model = SSIT;    % Create SSIT instance and call it 'Model'.

% Set species names for bursting gene expression model:
Model.species = {'x1';'x2';'x3'}; % Set species names for bursting gene expression model:
% x1 - OFF gene
% x2 - ON gene
% x3 - mRNA

% Set Stoichiometry of reactions:
Model.stoichiometry = [-1,1,0,0;...
                        1,-1,0,0;...
                        0,0,1,-1]; 

% Define a time-varying TF/MAPK input signal:
Model.inputExpressions = {'IHog','a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'};

% Set propensity functions:
Model.propensityFunctions = {'x1*IHog';'k21*x2';'kr*x2';'g*x3'}; 
% Model.propensityFunctions = {'k12*x1';'k21*x2*max(0,1-IHog)';'kr*x2';'g*x3'}; 

% Set initial condition:
Model.initialCondition = [1;0;0]; 

% Set times (s) at which to compute distributions:
Model.tSpan = [0:5:60]; 

%% Plot the TF/MAPK signal
% Next, we have to gues some initial guesses for parameters.
% First, let's tinker with the MAPK signal to get it to match somewhat
% qualitatively to what we see in experiments.  We don't have to be exact,
% ballpark parameters should be fine to start.
Model.parameters = ({'k21',30;'kr',100;'g',0.005; ...
    'a0',0.01;'a1',1;'r1',0.4;'r2',.1});
par = [Model.parameters{:,2}];
t = [0:60];
TF = par(4)+par(5)*exp(-par(6)*t).*(1-exp(-par(7)*t)).*(t>0);
figure(1); plot(t,TF,'linewidth',3); 
set(gca,'fontsize',16)
xlabel('Time (min)'); ylabel('Hog1(t)')
% Try tinkering with the MAPK signal parameters (parameters 5-8)to get it 
% to match somewhat qualitatively to what we see in experiments: maximum at
% ~2 minutes and adaptation in ~10 min.
% We don't have to be exact, ballpark parameters should be fine to start.

%% Solve and plot using the FSP approach
% To solve the model, we first select the solution scheme ('FSP') and then
% we call the SSIT.solve method.
Model.parameters = ({'k21',30;'kr',100;'g',0.005; ...
    'a0',0.01;'a1',1;'r1',0.4;'r2',.1});

Model.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
[FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis

% Next we make plots of the marginal distributions at time points 3, 5, 7,
% 9, 11, 13 and plot these in figures 1:3 for the three different species.
Model.makePlot(FSPsoln,'marginals',[3:2:13],false,(1:3))    % Plot marginal distributions

% We can also plot the means and standard deviations versus time in figure
% 100:
Model.makePlot(FSPsoln,'meansAndDevs',[],false,100)    % Plot marginal distributions

% Try to tune the parameters until you see:
% Bimodal expression (i.e., a population of active cells and a population of
% inactive cells).
% Perfect adaptation (all mRNA gone) at about 25 min.
% An average of ~50 mRNA at the highest expression time.

%% Load smFISH Data and compare to model
% In this section, we load some data to compare to the model.  For this
% example, we are going to use some data that Gregor Neuert collected.
Model = Model.loadData('../ExampleData/NeuertData/Result_Exp1_rep1_RNA_CY5_total_FORMATTED.csv',{'x3','mRNA'});

% The model should start at some unknown steady state distribution that
% depends on its parameters.  To reproduce this, we start the clock at a
% very early time (-120 min) to allow the model to relax to steady state
% before starting our analysis.
Model.initialTime = -120;
Model.tSpan = unique([Model.initialTime,Model.dataSet.times]);

% Next, we call the code to make the fitting plots.
Model.makeFitPlot

% After running this code, you will see a number of new plots:
% * Model and Data Means and Standard deviations versus time.
% * Model and Data probability mass functions versus time.
% * Model and Data cumulative distributions versus time.
% * Maximum Likelihood result versus time.
% The first three should be relatively self-explanatory.  The fourth one
% shows the likelihood function for each of the time points (blue) as well
% as the best possible likelihood if the model gave a perfect match (red
% line) and in orange what you might expect as a real fit for a perfectly
% identified model. In other words, a really good fit would be one where the blue line is
% close to the orange line, and the difference provides a sense as to how
% much the model could potentially be improved at each time point.

%% Fit the model to the smFISH data
% Once you have an okay guess for parameters, we can use this as an initial
% guess and let the computer try to identify better parameters.  Here, we
% will start by fitting on the first four parameters.
Model.fittingOptions.modelVarsToFit = [1:7];

% Here we use the current parameters as our initial guess:
x0 = [Model.parameters{Model.fittingOptions.modelVarsToFit,2}]';

% Here we call the search process with some fitting options.
fitOptions = optimset('Display','iter','MaxIter',1000);
[pars,likelihood] = Model.maximizeLikelihood(x0,fitOptions);

% Update Model and Make Plots of Results
Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(pars);
Model.makeFitPlot

% As the fit gets a little closer, you can also try to let the model fit
% the MAPK signal dynamics as well.  For the default data set,
% "Result_Exp1_rep1_RNA_CY5_total_FORMATTED", you should be able to get a
% fit with the MLE of better than 38500 after a few rounds of fitting. 

%% Does the Model predict a good MAPK(t)?
par = [Model.parameters{:,2}];
t = [0:60];
TF = par(4)+par(5)*exp(-par(6)*t).*(1-exp(-par(7)*t)).*(t>0);
figure(25); plot(t,TF,'linewidth',3); 
set(gca,'fontsize',16)
xlabel('Time (min)'); ylabel('Hog1(t)')


%% Quantifying model Sensitivities.
% By now, you have found a model that matches okay to your data.  (If not,
% you could add additional states or reactions to the system).  But just
% because you found one model that fits, does NOT mean that is the correct
% model.  There could be an infinite numbr of parameters that all match to
% the same data.  In this section, we are going to search around in
% parameter space to determine what is the uncertainty in the parameters
% given our model.

% In this first section, we are going to compute the sensitivity of the
% model to the different parameters.
Model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
Model.sensOptions.solutionMethod = 'finiteDifference';
[sensSoln] = Model.solve;  % Solve the sensitivity problem
Model.makePlot(sensSoln,'marginals',[],false,[11:13]) % Plot marginal sensitivities
% This will results in a few plots that show how chainging each of the
% model parameters would result in changes to the species' distributions.

% For later use, we are also going to compute the Fisher Information Matrix.  
fimResults = Model.computeFIM(sensSoln.sens);
[FIM] = Model.evaluateExperiment(fimResults,Model.dataSet.nCells);

% The inverse of the FIM provides an estimate of the model uncertainty.
% Here we are going to look at the FIM for the log of the model parameters
% and use that to compute the covariance of the log of the parameters.
% (Because parameters are positive values, but can very significantly in
% their magnitudes, it is often useful to examine them in a log-scale).
FIMlog = diag([Model.parameters{:,2}])*FIM*diag([Model.parameters{:,2}]);
covLog = FIMlog^-1;

%% Metropolis Hastings
% Now that we have an estimate of the shape of the uncertainty using the
% FIM we can now search parameter space and see what other parameter
% combinations are also closely matching to our data.  For this, we are
% going to use the Metropolis Hastings algorithm, where the proposal
% distribution is a multi-variate gaussian with a covariance that is
% proportional to the inverse FIM.  Here, we set up the MH parameters:
Model.solutionScheme = 'FSP'; % Set solutions scheme to FSP Sensitivity
Model.fittingOptions.modelVarsToFit = [1:7];
MHOptions = struct('numberOfSamples',1000,'burnin',0,'thin',3);
proposalWidthScale = 0.002;
MHOptions.proposalDistribution  = @(x)mvnrnd(x,proposalWidthScale*(covLog+covLog')/2);

% Next, we call the codes to sample the posterior parameter space:
[pars,likelihood,chainResults] = Model.maximizeLikelihood([],MHOptions,'MetropolisHastings');
% When this runs, you want to see an acceptance of about 0.3 to 0.4, meaning that
% about a third of the proposals are accepted.  If the number is too small
% you need to decrease the proposal width; if it is too large you may need
% to increase the proposal width. For the default data set and model, I
% found that a scale of .5 to 5% of the FIM-based COV led to an okay acceptance
% rate, but this is variable and will change depending on the initial value
% in the chain.

% Often the MH search can reveal a better parameter set, so let's make sure
% to update our model if it does:
Model.parameters(:,2) = num2cell(pars);
Model.makeFitPlot
% If you do notice better fits, it would be good to re-run the fminsearch
% again - it is possible you can still find a better model to explain your
% data.  This can take several rounds of iteration before convergence.  I
% recommend creating a while loop to make it automated.

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

% Next, we will show the scatter plot of a couple parameters.  It is
% helpful to show these in linear scale as well as in a natural log scale.
% For illustration, we also compare the spread of the posterior to the
% covariance predicted by the FIM from before.
Q = [3,4];
subplot(1,3,2)
Model.makeMleFimPlot(exp(chainResults.mhSamples)',FIM,Q,0.95,1); hold on
title('Posterior -- Linear Scale')
subplot(1,3,3)
Model.makeMleFimPlot(chainResults.mhSamples',FIMlog,Q,0.95,1); hold on
title('Posterior -- Log Scale')

%% Effective Sample Size
ipar = 5;
ac = xcorr(chainResults.mhSamples(:,ipar)-mean(chainResults.mhSamples(:,ipar)),'normalized');
ac = ac(size(chainResults.mhSamples,1):end);
plot(ac,'LineWidth',3)
N = size(chainResults.mhSamples,1);
tau = 1+2*sum(abs(ac(2:N/5)));
Neff = N/tau

%% Simulate Dataset
ModelSim = Model;
ModelSim.solutionScheme = 'SSA';  % Set solution scheme to SSA.
ModelSim.ssaOptions.Nexp = 1; 
ModelSim.ssaOptions.useTimeVar = true;
ModelSim.ssaOptions.signalUpdateRate = 5;
ModelSim.ssaOptions.nSimsPerExpt = 1000;
ModelSim.ssaOptions.applyPDO = false; % Include the distortion in the SSA data.
ModelSim.solve([],'HogSSAData.csv');   

%% Fit Simulated Data Set.
ModelSim.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
ModelSim = ModelSim.loadData('HogSSAData.csv',{'x3','exp1_s3'});
[FSPsoln,ModelSim.fspOptions.bounds] = ModelSim.solve;  % Solve the FSP analysis
ModelSim.makeFitPlot

%% Optimize an experiment
Model.parameters = ({'k21',30;'kr',100;'g',0.005; ...
    'a0',0.01;'a1',1;'r1',0.4;'r2',.1});
cellCounts = Model.dataSet.nCells;
nCellsTotal = sum(cellCounts);
% Compute FIM with all parameters free to vary
Model.fittingOptions.pdoVarsToFit = [];
Model.tSpan = [0:1:60]; 
Model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
Model.sensOptions.solutionMethod = 'finiteDifference';
[sensSoln] = Model.solve;  % Solve the sensitivity problem

fimResults = Model.computeFIM(sensSoln.sens); 
% FIMintuit = Model.evaluateExperiment(fimResults,cellCounts);
% Optimize cell counts for different objectives
nCellsOptDetA = Model.optimizeCellCounts(fimResults,nCellsTotal,'Determinant');
nCellsOptDetB = Model.optimizeCellCounts(fimResults,nCellsTotal,'[1:3]');
FIMOptA = Model.evaluateExperiment(fimResults,nCellsOptDetA);
FIMOptB = Model.evaluateExperiment(fimResults,nCellsOptDetB);
[diag(FIMOptA^-1),diag(FIMOptB^-1)]

