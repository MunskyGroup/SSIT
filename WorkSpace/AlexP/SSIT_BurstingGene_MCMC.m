clear all
close all
addpath(genpath('../src'));

%% Create Bursting Gene Model
F1 = SSIT;
F1.species = {'offGene';'onGene';'Protein'};
F1.initialCondition = [1;0;0];

F1.propensityFunctions = {'k1*offGene';'k2*onGene';...
    'k3*onGene';'k4*Protein'};
F1.stoichiometry = [-1,1,0,0;...
                         1,-1,0,0;...
                         0,0,1,-1];

F1.parameters = ({'k1',0.2;'k2',0.5;'k3',10;'k4',0.1}); 

%% Set priors. 
%% WARNING: 
%   Setting priors means that the "maximizeLikelihood" function 
%   will return the peak of the posterior, not the likelihood!  
%% Comment out the priors to return the max log-likelihood.

log10PriorMean = [-1 -0.5 0 -1];
log10PriorStd = 2*ones(1,4);

F1.fspOptions.initApproxSS = true;

F1.fittingOptions.modelVarsToFit = (1:4);
F1.fittingOptions.logPrior = @(x)-sum((log10(x)-log10PriorMean(1:4)).^2./(2*log10PriorStd(1:4).^2));


%% Solve using FSP
F1 = F1.formPropensitiesGeneral('BurstingGene');
F1.tSpan = [1:1:20];
F1.initialTime = 0;
F1.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
[FSPsoln,F1.fspOptions.bounds] = F1.solve;  % Solve the FSP analysis

%% Solve FSP Model again using the bounds from the last solution
% If we start with the bounds computed in the first analysis, 
% the solution is often much faster.
[FSPsoln] = F1.solve(FSPsoln.stateSpace);  % Solve the FSP analysis

%% Grab some data from FSP solution
F1.ssaOptions.nSimsPerExpt = 1; 
F1.ssaOptions.Nexp = 1;
% Set the seed (will always give same random data)
rng(123); 

dataFile = 'BurstingGene_data.csv';

F1.sampleDataFromFSP(FSPsoln,dataFile); % Sample some data from FSP solution

dataToFit = {'offGene','exp1_s1';'onGene','exp1_s2';'Protein','exp1_s3'};
%dataToFit = {'Protein','exp1_s3'};

% Add Data to Model
F1 = F1.loadData(dataFile,dataToFit);

fitParameters = [1:4];

% MLE Fitting Options
maxFitIter = 1000;
nFitRounds = 3;

% Metropolis Hastings Properties
nSamplesMH = 1000; % Number of MH Samples to run
nThinMH = 2; % Thin rate for MH sampling
nBurnMH = 100; % Number for MH burn in

% Number of MH samples to use for FIM calculation
nFIMsamples = 10;

% FIM options
fimScale = 'log'; % Maximize FIM for log parameters

% Plotting options
mhPlotScale = 'log10';  % Show MH and FIM plots in log10 scale.

%%
% Randomize the initial parameter set.
%np = size(F1.parameters,1);
%F1.parameters(:,2) = ...
%    num2cell([F1.parameters{:,2}]'.*(1+0.5*randn(np,1)));

% Fit Model to data.
fitOptions = optimset('Display','none','MaxIter',maxFitIter);
fitPars = F1.maximizeLikelihood([],fitOptions);
fitPars

%% Compute FIM
%F2 = F1.computeFIM;

%% Run MH on GR Models.
MHFitOptions.thin=nThinMH;
MHFitOptions.numberOfSamples=nSamplesMH;
MHFitOptions.burnIn=nBurnMH;
MHFitOptions.progress=true;
MHFitOptions.numChains = 1;
%MHFitOptions.useFIMforMetHast = true;
MHFitOptions.logForm = false;
MHFitOptions.proposalDistribution=@(x)abs(x+0.1*randn(size(x)));
%F1.fittingOptions.modelVarsToFit = fitParameters;
MHFitOptions.saveFile = 'BurstingGene_MHtrace.mat';
delete 'BurstingGene_MHtrace.mat'
[~,~,MHResults] = F1.maximizeLikelihood(...
    [], MHFitOptions, 'MetropolisHastings');
%[~,~,MHResults] = F2.maximizeLikelihood(...
    %fitPars, MHFitOptions, 'MetropolisHastings');
%delete(MHFitOptions.saveFile)