%% dataFittingExample
% In this script, we show how to load a data set and then fit it using
% the FSP tools.
clear all
addpath(genpath('../src'));
%% Load a model template
ModelChoice = 'BurstingGene';  % Two species problem (mRNa and protein)
F2 = SSIT(ModelChoice);

%% Adjust to add time-varying reaction with new parameters.
F2.propensityFunctions{1} = 'kon*I1*(2-x1)';
F2.inputExpressions = {'I1','a0+a1*exp(-r1*t)*(1-exp(-r2*t))^eta*(t>0)'};
F2.parameters = ({'koff',0.14; ...
    'kon',0.14; ...
    'kr',25; ...
    'gr',0.01; ...
    'a0',0.006; ...
    'a1',0.4; ...
    'r1',0.04; ...
    'r2',0.1; ...
    'eta',1});

%% Load some data
linkedSpecies = {'x2','RNA_nuc'}; % link species x2 to the column 'RNA_nuc'
F2 = F2.loadData('../examples/Example/DUSP1_Dex_100nM_Rep1_Rep2.csv',linkedSpecies);
F2.initialTime = -120;  % set initial time to -120 (simulate approach to equilibrium before input)
F2.tSpan = unique([F2.initialTime,F2.dataSet.times]);

%% Compute Likelihood of Data given the FSP model 
F2.solutionScheme = 'FSP';    % Set solution scheme to FSP.
[FSPsoln,bounds] = F2.solve;  % Solve the FSP analysis
F2.fspOptions.bounds = bounds;% Save bound for faster analyses 
F2.pdoOptions.PDO = []; % Assume perfect measurements (see below for relaxed example).
%%
pars = [F2.parameters{:,2}]';
[fit_error,fitSolutuions] = F2.computeLikelihood(pars,FSPsoln.stateSpace);

%% Make plot of fitting results
F2.makeFitPlot(fitSolutuions)

%% Run a search algorithm to maximize the likelihood function
pars = [F2.parameters{:,2}]';
x0 = log10(pars);
obj = @(x)-F2.computeLikelihood(10.^x,FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.
fitOptions = optimset('Display','iter','MaxIter',20);
obj(x0)
for i=1:3
    x0  = fminsearch(obj,x0,fitOptions);
    F2.parameters(:,2) = num2cell(10.^x0);
    [FSPsoln,F2.fspOptions.bounds] = F2.solve;  % Solve the FSP analysis
    obj = @(x)-F2.computeLikelihood(10.^x,FSPsoln.stateSpace);  % Update the statespace for next round.
end
[fit_error,fitSolutuions] = F2.computeLikelihood(pars,FSPsoln.stateSpace);
F2.makeFitPlot(fitSolutuions)

%% Fit model to data given given a prior on parameters
% Here we will implement a lognormal prior on parameters where all
% pararameters are expected to be 1 with a log10 deviation of 1.
mu = zeros(9,1); 
mu(4) = -2; % degradation rate is probably on order of 10^-2
sig = ones(9,1);
F3 = F2;
F3.fittingOptions.logPrior = @(p)-(log10(p)-mu).^2./(2*sig);
obj = @(x)-F3.computeLikelihood(10.^x,FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.
fitOptions = optimset('Display','iter','MaxIter',20);
x0 = log10([F3.parameters{:,2}]');
for i=1:3
    x0  = fminsearch(obj,x0,fitOptions);
    F3.parameters(:,2) = num2cell(10.^x0);
    [FSPsoln,F3.fspOptions.bounds] = F3.solve;  % Solve the FSP analysis
    obj = @(x)-F3.computeLikelihood(10.^x,FSPsoln.stateSpace);  % Update the statespace for next round.
end
[fit_error,fitSolutuions] = F3.computeLikelihood(pars,FSPsoln.stateSpace);
F3.makeFitPlot(fitSolutuions)

%% Fit model to data given the FSP model and a parametrized PDO
F4=F3;
F4.solutionScheme = 'FSP';    % Set solution scheme to FSP.
[FSPsoln,bounds] = F4.solve;  % Solve the FSP analysis
F4.fspOptions.bounds = bounds;% Save bound for faster analyses 

mu = zeros(9,1); 
mu(4) = -2; % degradation rate is probably on order of 10^-2
sig = ones(9,1);
F4.fittingOptions.logPrior = @(p)[-(log10(p(1:9))-mu).^2./(2*sig);...
    -max(0,1000*(p(10)-1))];

% Define probabilistic distortion for each species S1, S2,...
F4.pdoOptions.type = 'Binomial - Parametrized';
F4.pdoOptions.props.CaptureProbabilityS1 = @(x,C)0;  % Use zero for unobserved species.
F4.pdoOptions.props.CaptureProbabilityS2 = @(x,C)1 - C(1)*x/(40+x)'; % Use 1.0 if there is no loss.
F4.pdoOptions.props.ParameterGuess = 0.000001; 
F4.fittingOptions.pdoVarsToFit = 'all';

pars = [[F4.parameters{:,2}]';F4.pdoOptions.props.ParameterGuess];
[fit_error] = F4.computeLikelihood(pars,FSPsoln.stateSpace)
