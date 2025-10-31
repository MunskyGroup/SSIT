close all
clear
addpath(genpath('../src'));

rng(1)  % Set RNG seed for reproducibility.
%% Define Model
% example = 'Poisson';
example = 'Poisson';
switch example
    case 'Poisson'
        ModelTrue = SSIT;
        ModelTrue.species = {'rna'};
        ModelTrue.initialCondition = 0;
        ModelTrue.propensityFunctions = {'kr';'gr*rna'};
        ModelTrue.stoichiometry = [1,-1];
        ModelTrue.fspOptions.fspTol = 1e-5;
        ModelTrue.parameters = ({'kr',10;'gr',0.3});
        ModelTrue = ModelTrue.formPropensitiesGeneral('Poiss',true);
        dataToFit = {'rna','exp1_s1'};
        fitParameters = [1:2];
    case 'Toggle'
        ModelTrue = SSIT;
        ModelTrue.parameters = {'kb',10;'ka',80;'M',20;'g',1};
        ModelTrue.species = {'lacI';'lambdaCI'};
        ModelTrue.stoichiometry = [1,-1,0, 0;...
            0, 0,1,-1];
        ModelTrue.propensityFunctions = {'kb+ka*M^3/(M^3+lambdaCI^3)';...
            'g*lacI';...
            'kb+ka*M^3/(M^3+lacI^3)';...
            'g*lambdaCI'};
        ModelTrue.initialCondition = [0;0];
        ModelTrue.customConstraintFuns = {'(x1-3).^2.*(x2-3).^2'};
        dataToFit = {'lacI','exp1_s1';'lambdaCI','exp1_s2'};
        fitParameters = [1:4];
end

saveFileName = ['IterativeExperimentResults_',example];

%% Generate Model Propensity Functions and Solve True Model
ModelTrue = ModelTrue.formPropensitiesGeneral('TrueModel',true);
[ModelSolution,ModelTrue.fspOptions.bounds] = ModelTrue.solve;

%% Set Rules for Experiment Designs and Iterations
% Possible Time points for Experiments.
nT = 21;
ModelTrue.tSpan = linspace(0,10,nT);

% Number of experiment Rounds
nExptRounds = 5;
    
% MLE Fitting Options
maxFitIter = 1000;
nFitRounds = 3;

% Metropolis Hastings Properties
nSamplesMH = 1000; % Number of MH Samples to run
nThinMH = 2; % Thin rate for MH sampling
nBurnMH = 100; % Numbr for MH burn in

% Number of MH samples to use for FIM calculation
nFIMsamples = 10;

% Total number of new cells allowed in each experiment
numCellsPerExperiment = 30; 

% Definition of initial experiment
nextExperiment = zeros(1,nT);
nextExperiment([1,11,21]) = round(numCellsPerExperiment/3);

% FIM options
fimScale = 'log'; % Maximize fim for log parameters

% Plotting options
mhPlotScale = 'log10';  % Show MH and FIM plots in log10 scale.

%%
% Create Model for Estimate.
ModelGuess = ModelTrue;

% Randomize the initial parameter set.
np = size(ModelGuess.parameters,1);
ModelGuess.parameters(:,2) = ...
    num2cell([ModelGuess.parameters{:,2}]'.*(1+0.5*randn(np,1)));

% Initialize cells for saved results
covMH = cell(1,nExptRounds);
covLogMH = cell(1,nExptRounds);
covLogFIM_Prediction = cell(1,nExptRounds);
covFIM_Prediction = cell(1,nExptRounds);
parametersFound = cell(1,nExptRounds);
exptDesigns = cell(1,nExptRounds);
FIMOptNextExptSaved = cell(1,nExptRounds);
MHResultsSaved = cell(1,nExptRounds);

% Loop over subsequent experiments
allDataSoFar = [];
nTotalCells = zeros(1,nT);
for iExpt = 1:nExptRounds
    % Generate "Fake" data
    dataFile = ['FakeExperiment_',num2str(iExpt),'.csv'];
    ModelTrue.ssaOptions.nSimsPerExpt = max(nextExperiment);
    ModelTrue.ssaOptions.Nexp = 1;
    ModelTrue.sampleDataFromFSP(ModelSolution,dataFile)
    
    % DownSelect fake data to specified number of cells at each time
    X = importdata(dataFile);
    for iT = 1:length(ModelTrue.tSpan)
        if nextExperiment(iT)>0
            J = find(X.data(:,1)==ModelTrue.tSpan(iT),nextExperiment(iT));
            
            % Add new data to previous data set
            allDataSoFar = [allDataSoFar;X.data(J,:)];
        end
    end

    % Save generated data to file.
    A = table;
    for j=1:length(X.colheaders)
        A.(X.colheaders{j}) = allDataSoFar(:,j);
        writetable(A,dataFile)
    end
    
    nTotalCells = nTotalCells + nextExperiment

    % Add Data to Model
    ModelGuess = ModelGuess.loadData(dataFile,dataToFit);
            
    % Fit Model to data.
    fitOptions = optimset('Display','none','MaxIter',maxFitIter);
    for iFitIter=1:nFitRounds
        fitPars = ModelGuess.maximizeLikelihood([],fitOptions);
        ModelGuess.parameters(:,2) = num2cell(fitPars);
    end

    % Run Metropolis Hastings
    MHFitOptions.thin=nThinMH;
    MHFitOptions.numberOfSamples=nSamplesMH;
    MHFitOptions.burnIn=nBurnMH;
    MHFitOptions.progress=true;
    MHFitOptions.useFIMforMetHast=true;
    MHFitOptions.CovFIMscale = 1.0;
    MHFitOptions.numChains = 1;
    MHFitOptions.saveFile = 'TMPMHChain.mat';
    delete 'TMPMHChain.mat'
    ModelGuess.fittingOptions.modelVarsToFit = fitParameters;
    [newPars,~,MHResults] = ModelGuess.maximizeLikelihood(...
        [], MHFitOptions, 'MetropolisHastings');
    ModelGuess.parameters(ModelGuess.fittingOptions.modelVarsToFit,2) = num2cell(newPars);
    delete('TMPMHChain.mat')

    % Compute FIM for subsampling of MH results.
    J = floor(linspace(nSamplesMH/2,nSamplesMH,nFIMsamples));
    MHSamplesForFIM = exp(MHResults.mhSamples(J,:));
    ModelGuess.tSpan = ModelTrue.tSpan;
    fimResults = ModelGuess.computeFIM([],fimScale,MHSamplesForFIM);

    % Find optimal NEXT experiment design given parameter sets
    nextExperiment = ModelGuess.optimizeCellCounts(fimResults,numCellsPerExperiment,'Determinant',[],nTotalCells);
    FIMOptNextExpt = ModelGuess.totalFim(fimResults,nextExperiment+nTotalCells);

    % Compute and Save Covariance from MH from current stage.
    parametersFound{iExpt} = newPars;
    FIMOptNextExptSaved{iExpt} = FIMOptNextExpt;
    covLogMH{iExpt} = cov(MHResults.mhSamples);
    covMH{iExpt} = cov(exp(MHResults.mhSamples));
    MHResultsSaved{iExpt} = MHResults;

    % Save Predicted FIM for next stage.
    for jFIM=1:nFIMsamples
        covFIM_Prediction{iExpt}{jFIM} = inv(FIMOptNextExpt{jFIM});
    end

    exptDesigns{iExpt} = nextExperiment;

    f = figure;
    f.Name = ['Current MH Results and Next FIM Prediction (Round ',num2str(iExpt),')'];
    ModelGuess.plotMHResults(MHResults,FIMOptNextExpt,fimScale,mhPlotScale,f)
    if iExpt>1
        f = figure;
        f.Name = ['Current MH Results and Previous FIM Prediction (Round ',num2str(iExpt),')'];
        ModelGuess.plotMHResults(MHResults,[FIMOptNextExptSaved{iExpt-1}],fimScale,mhPlotScale,f)
    end

    save(saveFileName,'parametersFound','FIMOptNextExptSaved','covMH',...
        'covLogMH','exptDesigns','')

end


