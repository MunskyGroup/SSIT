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

F1.parameters = ({'k1',0.2;'k2',0.5;'k3',10;'k4',0.1}); % TODO: check!

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

%% WARNING:
%   Will operating on log10(x) parameters mess up the Hastings ratio?  
%   Will it be biased for large x?

%% (2) Solving a Model
%%      (2A) Solve using FSP
F1 = F1.formPropensitiesGeneral('BurstingGene');
F1.tSpan = [1:1:20];
F1.initialTime = 0;
F1.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
[FSPsoln,F1.fspOptions.bounds] = F1.solve;  % Solve the FSP analysis

%%          (2A.1) Solve FSP Model again using the bounds from the last solution
% If we start with the bounds computed in the first analysis, the solution is
% often much faster.
[FSPsoln] = F1.solve(FSPsoln.stateSpace);  % Solve the FSP analysis

%% Grab some data from FSP solution

F1.ssaOptions.nSimsPerExpt = 1; 
F1.ssaOptions.Nexp = 1;
rng(123); % Set the seed (will always give same random data)

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
nBurnMH = 100; % Numbr for MH burn in

% Number of MH samples to use for FIM calculation
nFIMsamples = 10;

% FIM options
fimScale = 'log'; % Maximize fim for log parameters

% Plotting options
mhPlotScale = 'log10';  % Show MH and FIM plots in log10 scale.

%%
% Randomize the initial parameter set.
%np = size(F1.parameters,1);
%F1.parameters(:,2) = ...
%    num2cell([F1.parameters{:,2}]'.*(1+0.5*randn(np,1)));

% Fit Model to data.
fitOptions = optimset('Display','none','MaxIter',maxFitIter);
%for iFitIter=1:nFitRounds
      fitPars = F1.maximizeLikelihood([],fitOptions);
%     F1.parameters(:,2) = num2cell(fitPars);
%end

fitPars

%% Compute FIM
%F2 = F1.computeFIM;

%%    Run MH on GR Models.
MHFitOptions.thin=nThinMH;
MHFitOptions.numberOfSamples=nSamplesMH;
MHFitOptions.burnIn=nBurnMH;
MHFitOptions.progress=true;
MHFitOptions.numChains = 1;
%MHFitOptions.useFIMforMetHast = true;
%MHFitOptions.logForm = false;
%F1.fittingOptions.modelVarsToFit = fitParameters;
MHFitOptions.saveFile = 'BurstingGene_MHtrace.mat';
delete 'BurstingGene_MHtrace.mat'
[~,~,MHResults] = F1.maximizeLikelihood(...
    [], MHFitOptions, 'MetropolisHastings');
%[~,~,MHResults] = F2.maximizeLikelihood(...
    %fitPars, MHFitOptions, 'MetropolisHastings');
%delete(MHFitOptions.saveFile)

% Compute FIM for subsampling of MH results.
    %J = floor(linspace(nSamplesMH/2,nSamplesMH,nFIMsamples));
    %MHSamplesForFIM = exp(MHResults.mhSamples(J,:));
    %F1.tSpan = ModelTrue.tSpan;
    %fimResults = F1.computeFIM([],fimScale,MHSamplesForFIM);

    %f = figure;
    %f.Name = ['Current MH Results and Next FIM Prediction (Round ',num2str(iExpt),')'];
    %F1.plotMHResults(MHResults,FIMOptNextExpt,fimScale,mhPlotScale,f)
    %if iExpt>1
    %    f = figure;
    %    f.Name = ['Current MH Results and Previous FIM Prediction (Round ',num2str(iExpt),')'];
    %    F1.plotMHResults(MHResults,[FIMOptNextExptSaved{iExpt-1}],fimScale,mhPlotScale,f)
    %end

    %save(saveFileName,'parametersFound','FIMOptNextExptSaved','covMH',...
    %    'covLogMH','exptDesigns','')

%%          (2A.2) Make plots of FSP solution
%F1.makePlot(FSPsoln,'meansAndDevs',[],[],1,{'linewidth',3,'color',[0,1,1]}) % Make plot of mean vs. time.
%F1.makePlot(FSPsoln,'marginals',[],[],2,{'linewidth',3,'color',[0,0,1]}) % Make plot of mean vs. time.

%%      (2B) Solve using SSA
%F2 = F1;
%F2.solutionScheme = 'SSA';
%SSASoln = F2.solve;

%%          (2B.1) Make plots of SSA solution
%F2.makePlot(SSASoln,'trajectories',[],[],4) % Make some plots.
%F1.makePlot(FSPsoln,'meansAndDevs',[],[],4,...
%    {'linewidth',4,'color',[0,1,1],'Marker','s','MarkerSize',20}) % Add FSP Solution to plot.

