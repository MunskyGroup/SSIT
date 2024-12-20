addpath(genpath('../../src/Commandline/'));

function F1 = setupAndSolveModel()
    %% Create Bursting Gene Model
    F1 = SSIT;
    F1.species = {'offGene'; 'onGene'; 'Protein'};
    F1.initialCondition = [1; 0; 0];

    F1.propensityFunctions = {'k1*offGene'; 'k2*onGene'; 'k3*onGene'; 'k4*Protein'};
    F1.stoichiometry = [-1,1,0,0; 1,-1,0,0; 0,0,1,-1];
    F1.parameters = ({'k1',0.2; 'k2',0.5; 'k3',10; 'k4',0.1}); 
    F1.fspOptions.initApproxSS = true;

    %% Solve using FSP
    F1 = F1.formPropensitiesGeneral('BurstingGene');
    F1.tSpan = [1:0.1:5];
    F1.initialTime = 0;
    F1.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
    [FSPsoln, F1.fspOptions.bounds] = F1.solve;  % Solve the FSP analysis

    %% Solve FSP Model again using the bounds from the last solution
    % If we start with the bounds computed in the first analysis, 
    % the solution is often much faster.
    [FSPsoln] = F1.solve(FSPsoln.stateSpace);  % Solve the FSP analysis

    %% Grab some data from FSP solution
    F1.ssaOptions.nSimsPerExpt = 1; 
    F1.ssaOptions.Nexp = 1;
    % Set the seed (will always give same random data)
    rng(123); 

    % Sample some data from FSP solution
    %dataFile = 'BurstingGene_data_FSP.csv';
    %F1.sampleDataFromFSP(FSPsoln, dataFile); 

    %% Alt: Read in SSA sample used in Stan (for comparison)
    dataFile = 'BurstingGene_data_SSA.csv';

    dataToFit = {'offGene','exp1_s1';'onGene','exp1_s2';'Protein','exp1_s3'};

    % Add Data to Model
    F1 = F1.loadData(dataFile, dataToFit);
end

function result = performMLEandMCMC(F1,logSpace,stepSize,chainLength)

    % MLE Fitting Options
    maxFitIter = 1000;
    nFitRounds = 3;

    %% Set priors. 
    %% WARNING: 
    %   Setting priors means that the "maximizeLikelihood" function 
    %   will return the peak of the posterior, not the likelihood!  
    %% Comment out the priors to return the max log-likelihood.
    log10PriorMean = [-1 -0.5 1 -1];
    log10PriorStd = 2*ones(1,4);

    F1.fittingOptions.modelVarsToFit = (1:4);
    F1.fittingOptions.logPrior = @(x)-sum((log10(x)-log10PriorMean(1:4)).^2./(2*log10PriorStd(1:4).^2));

    % Metropolis Hastings Properties
    nSamplesMH = chainLength; % Number of MH Samples to run
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
    % Custom initialization using log-normal distribution
    np = size(F1.parameters, 1);
    if logSpace
        F1.parameters(:, 2) = num2cell(exp(randn(np, 1)));  % log-normal(0, 1)
    else
        F1.parameters(:,2) = ...
        num2cell([F1.parameters{:,2}]'.*(abs(randn(np,1)))); % normal - CHECK
    end

    % Fit Model to data.
    fitOptions = optimset('Display','none','MaxIter',maxFitIter);
    fitPars = F1.maximizeLikelihood([],fitOptions);
    disp('Fitted Parameters:');
    disp(fitPars);

    %% Compute FIM
    % F2 = F1.computeFIM;

    %% Run MH on GR Models.
    MHFitOptions.thin = nThinMH;
    MHFitOptions.numberOfSamples = nSamplesMH;
    MHFitOptions.burnIn = nBurnMH;
    MHFitOptions.progress = true;
    MHFitOptions.numChains = 1;
    % MHFitOptions.useFIMforMetHast = true;

    % logForm defaults to true; 
    % for operating on parameters in linear space, use absolute values to
    % reflect negative proposals generated by randn across 0 boundary
    MHFitOptions.logForm = logSpace;  
    if ~logSpace
        MHFitOptions.proposalDistribution = @(x)abs(x+stepSize*randn(size(x)));
    end
    MHFitOptions.saveFile = 'BurstingGene_MHtrace.mat';
    delete 'BurstingGene_MHtrace.mat';
    [~,~,MHResults] = F1.maximizeLikelihood([], MHFitOptions, 'MetropolisHastings');
    %[~,~,MHResults] = F2.maximizeLikelihood(...
    %fitPars, MHFitOptions, 'MetropolisHastings');
    %delete(MHFitOptions.saveFile)
    disp('MH Results:');
    result = MHResults;
end

%% ChatGPT Plot MH Traces
function plotMHTraces(samples, paramNames)
    % Check if the number of parameter names matches the number of columns in samples
    if length(paramNames) ~= size(samples, 2)
        error('The number of parameter names must match the number of columns in samples.');
    end

    % Number of parameters
    numParams = size(samples, 2);

    % Create a figure for the trace plots
    figure;

    % Plot each parameter's trace
    for i = 1:numParams
        subplot(numParams, 1, i); % Create a subplot for each parameter
        plot(samples(:, i));      % Plot the trace of the parameter
        xlabel('Iteration');
        ylabel(paramNames{i});
        title(['Trace plot for ', paramNames{i}]);
    end

    % Set a global title for all subplots
    subplot(numParams, 1, 1);
    sgtitle('Metropolis-Hastings Trace Plots');
end

function ess = effective_sample_size(mhResults)
    ac = xcorr(mhResults.mhValue-mean(mhResults.mhValue),'normalized');
    ac = ac(size(mhResults.mhValue,1):end);
    % Uncomment following plot line to plot autocorrelation function,
    %% BUT this will overwrite one of the trace plots.
    % plot(ac,'LineWidth',3); hold on
    N = size(mhResults.mhValue,1);
    tau = 1+2*sum((ac(2:N/100)));
    Neff = N/tau;
    ess = Neff;
end

F1 = setupAndSolveModel();

%% Get computation time for MCMC
% linear space - careful! change the intialization!
%tic
%MHResults = performMLEandMCMC(F1,false,0.1,15000);
%toc

% log space
tic
logMHResults = performMLEandMCMC(F1,true,0.1,15000);
toc

transformed_logMHSamples = exp(logMHResults.mhSamples);

% Assume MHResults.mhSamples is a 1000x4 matrix
samples = MHResults.mhSamples;

% Define parameter names
paramNames = {'offGene', 'onGene', 'translation', 'degradation'};

% Call the function to plot the traces
plotMHTraces(samples, paramNames);

ess = effective_sample_size(MHResults);
ess

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