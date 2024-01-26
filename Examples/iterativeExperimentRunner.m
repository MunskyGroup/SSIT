function iterativeExperimentRunner(example,data,sampleType,nExptRounds,rngSeed,ind)
arguments
    example = 'Poisson'
    data = 'simulated';
    sampleType = 'random';
    nExptRounds = 5;
    rngSeed = []
    ind = randi(10000)
end

addpath(genpath('../src'));
if ~isempty(rngSeed)
    rng(rngSeed)
end

%% Define Model
% Number of experiment Rounds

showPlots = true;

switch example
    case 'Poisson'
        ModelTrue = SSIT;
        ModelTrue.species = {'rna'};
        ModelTrue.initialCondition = 0;
        ModelTrue.propensityFunctions = {'kr';'gr*rna'};
        ModelTrue.stoichiometry = [1,-1];
        ModelTrue.fspOptions.fspTol = 1e-5;
        ModelTrue.parameters = ({'kr',10;'gr',0.3});
        ModelTrue = ModelTrue.formPropensitiesGeneral(['Poiss_',num2str(ind)],true);
        dataToFit = {'rna','exp1_s1'};
        fitParameters = [1:2];
        nT = 21;
        ModelTrue.tSpan = linspace(0,10,nT);
        fimMetric = 'Determinant';
        nSamplesMH = 2000; % Number of MH Samples to run
        
        % Total number of new cells allowed in each experiment
        numCellsPerExperiment = 50;

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
        nT = 21;
        ModelTrue.tSpan = linspace(0,10,nT);
        fimMetric = 'Determinant';
        
        nSamplesMH = 2000; % Number of MH Samples to run

        % Total number of new cells allowed in each experiment
        numCellsPerExperiment = 50;
        ModelTrue = ModelTrue.formPropensitiesGeneral(['Toggle_',num2str(ind)],true);

    case 'DUSP1'
        ModelTrue = SSIT; % Rename SSIT code to make changes with out changing the original code
        ModelTrue.species = {'x1';'x2'}; % x1: number of on alleles,  x2: mRNA 
        ModelTrue.initialCondition = [0;0]; % Set initial conditions
        ModelTrue.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'}; % Set the propensity functions of the 4 reactions
        % Associated stoichiometry of the four reactions
        stoich = [1,0;   % Gene allele switching to on state
                 -1,0;   % Gene allele switching to off state
                  0,1;   % mRNA production
                  0,-1]; % mRNA degredation   
        ModelTrue.stoichiometry = transpose(stoich);
        % Input expression for time varying signals
        ModelTrue.inputExpressions = {'IGR','kcn0/knc+(t>=0)*kcn1/(r1-knc)*(exp(-knc*t)-exp(-r1*t))'};       
        % Defining the values of each parameter used
        ModelTrue.parameters = ({'koff',0.14;'kon',1;'kr',1;'gr',0.01;...
                           'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
%         ModelTrue.parameters=load('SGRS_model_v1.mat').SGRS_Model.parameters; % True model parameters after fitting to data

        ModelTrue.fspOptions.initApproxSS = true;
        
        dataToFit = {'x2','exp1_s2'};
        fitParameters = [1:4];
        timeDUSP1 = [0 10 20 30 40 50 60 75 90 120 150 180];
        nT = length(timeDUSP1);
        ModelTrue.tSpan = timeDUSP1;
        ModelTrue.fspOptions.bounds(3:4) = [2.1,100];
        fimMetric = 'TR1:4';
        ModelTrue.pdoOptions.unobservedSpecies = {'x1'};

        muLog10Prior = [-1 -1 -1 -2 -2 -2 -2 -2];
        sigLog10Prior = [2 2 2 2 2 2 2 2];
        ModelTrue.fittingOptions.modelVarsToFit = length(sigLog10Prior);
        muLog10Prior = muLog10Prior(ModelTrue.fittingOptions.modelVarsToFit);
        sigLog10Prior = sigLog10Prior(ModelTrue.fittingOptions.modelVarsToFit);
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);

        nSamplesMH = 2000; % Number of MH Samples to run

%         intuitiveDesign = [50 0 50 0 50 0 50 0 50 0 0 50;
%                              75 0 0 75 0 0 75 0 0 0 75 0;
%                              0 0 0 100 0 0 100 0 0 100 0 0;
%                              0 0 0 100 0 0 0 100 0 100 0 0;
%                              0 0 0 0 100 0 100 0 0 100 0 0];

                    intuitiveDesign = [20 0 20 0 20 0 20 0 20 0 0 20;
                                     30 0 0 30 0 0 30 0 0 0 30 0;
                                     0 0 0 40 0 0 40 0 0 40 0 0;
                                     0 0 0 40 0 0 0 40 0 40 0 0;
                                     0 0 0 0 40 0 40 0 0 40 0 0];

        % Total number of new cells allowed in each experiment
        numCellsPerExperiment = 120;
        ModelTrue = ModelTrue.formPropensitiesGeneral(['DUSP1_',num2str(ind)],true);

end
ModelTrue.fspOptions.fspTol = 1e-6;
saveFileName = ['IterativeExperimentResults_',example];

%% Generate Model Propensity Functions and Solve True Model
% ModelTrue = ModelTrue.formPropensitiesGeneral('TrueModel',true);
ModelTrue.fspOptions.verbose = true;
[ModelSolution,ModelTrue.fspOptions.bounds] = ModelTrue.solve;

%% FIM options
fimScale = 'log'; % Maximize fim for log parameters

%% True Model FIM
fimTrue = ModelTrue.computeFIM([],fimScale);

%% Verify that the true model and simulated data look correct.
if showPlots
    dataFile = ['FakeExperiment_TMP_',num2str(ind),'.csv'];
    nextExperiment = 100*ones(1,nT)
    ModelTrue.ssaOptions.nSimsPerExpt = max(nextExperiment);
    ModelTrue.ssaOptions.Nexp = 1;
    ModelTrue.sampleDataFromFSP(ModelSolution,dataFile)
    ModelTMP = ModelTrue.loadData(dataFile,dataToFit);
    ModelTMP.makeFitPlot
end

%% Set Rules for Experiment Designs and Iterations
% Possible Time points for Experiments.
    
% MLE Fitting Options
maxFitIter = 2000;
nFitRounds = 3;

% Metropolis Hastings Properties
nThinMH = 2; % Thin rate for MH sampling
nBurnMH = 100; % Number for MH burn in

% Number of MH samples to use for FIM calculation
nFIMsamples = 10;

% Definition of initial experiment
nextExperiment = zeros(1,nT);
nextExperiment([1,round(nT/2),nT]) = round(numCellsPerExperiment/3);

% Plotting options
mhPlotScale = 'log10';  % Show MH and FIM plots in log10 scale.

%% Experiment Design Definitions
incrementAdd = 5; % group size for Cells added
N = numCellsPerExperiment/incrementAdd; % needs to be int
randomCell = zeros(nExptRounds,nT);
for i = 1:nExptRounds
    n = randi(nT,1,N);
    for j = n
        randomCell(i,j) = randomCell(i,j) + incrementAdd;
    end
end

uniformCell = floor(ones(nExptRounds,nT)*numCellsPerExperiment/nT);
for i=1:nExptRounds
    if sum(uniformCell(i,:))<numCellsPerExperiment
        uniformCell(i,1:numCellsPerExperiment-sum(uniformCell(i,:)))=...
            uniformCell(i,1:numCellsPerExperiment-sum(uniformCell(i,:)))+1;
    end
end
%%
% Create Model for Estimate.
ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
ModelGuess = ModelTrue;
% ModelGuess.fspOptions.fspTol = inf;

% Randomize the initial parameter set.
np = length(fitParameters);
ModelGuess.parameters(fitParameters,2) = ...
    num2cell([ModelGuess.parameters{fitParameters,2}]'.*(1+0.5*randn(np,1)));

% Initialize cells for saved results
covMH = cell(1,nExptRounds);
covLogMH = cell(1,nExptRounds);
covLogFIM_Prediction = cell(1,nExptRounds);
covFIM_Prediction = cell(1,nExptRounds);
parametersFound = cell(1,nExptRounds);
exptDesigns = cell(1,nExptRounds);
% FIMOptNextExptSaved = cell(1,nExptRounds);
FIMcurrentExptSaved = cell(1,nExptRounds);
FIMcurrentExptTrueSaved = cell(1,nExptRounds);

MHResultsSaved = cell(1,nExptRounds);

% Loop over subsequent experiments
allDataSoFar = [];
nTotalCells = zeros(1,nT);  

for iExpt = 1:nExptRounds
    

    nTotalCells = nTotalCells + nextExperiment
    
    switch data
        case 'real'
            % Sample from real data
            [simData,csvFile] = sampleExperimentSim('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',timeDUSP1,nTotalCells);
            dataFile = csvFile;
            dataToFit = {'x2','RNA_nuc'};
        case 'simulated'
            % Generate "Fake" data
            dataFile = ['FakeExperiment_',num2str(iExpt),'_',num2str(ind),'.csv'];
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
    end

    % Add Data to Model
    ModelGuess = ModelGuess.loadData(dataFile,dataToFit);
            
    % Fit Model to data.
    fitOptions = optimset('Display','none','MaxIter',maxFitIter);
    for iFitIter=1:nFitRounds
        fitPars = ModelGuess.maximizeLikelihood([],fitOptions);
        ModelGuess.parameters(fitParameters,2) = num2cell(fitPars);
    end
    
    if showPlots
        ModelGuess.makeFitPlot
    end

    % Run Metropolis Hastings
    MHFitOptions.thin=nThinMH;
    MHFitOptions.numberOfSamples=nSamplesMH;
    MHFitOptions.burnIn=nBurnMH;
    MHFitOptions.progress=true;
    MHFitOptions.useFIMforMetHast=true;
    MHFitOptions.CovFIMscale = 0.95;
    MHFitOptions.numChains = 1;
    MHFitOptions.saveFile = ['TMPMHChain_',num2str(ind),'.mat'];
    delete(MHFitOptions.saveFile) 
    [newPars,~,MHResults] = ModelGuess.maximizeLikelihood(...
        [], MHFitOptions, 'MetropolisHastings');
    ModelGuess.parameters(ModelGuess.fittingOptions.modelVarsToFit,2) = num2cell(newPars);
    delete(MHFitOptions.saveFile)

    % Compute FIM for subsampling of MH results.
    J = floor(linspace(nSamplesMH/2,nSamplesMH,nFIMsamples));
    MHSamplesForFIM = exp(MHResults.mhSamples(J,:));
    ModelGuess.tSpan = ModelTrue.tSpan;

%     ModelGuess.evaluateExperiment(fimResults,nTotalCells,SIGprior); %%%%%%%%%%%%%%%%%%%
    fimResults = ModelGuess.computeFIM([],fimScale,MHSamplesForFIM);

    % FIM current experiment
    FIMCurrentExpt = ModelGuess.totalFim(fimResults,nTotalCells);

    % True FIM for chosen experiment
    FIMCurrentExpt_True = ModelTrue.totalFim(fimTrue,nTotalCells);
    
    switch sampleType
        case 'FIMopt'
            % Find optimal NEXT experiment design given parameter sets
            nextExperiment = ModelGuess.optimizeCellCounts(fimResults,numCellsPerExperiment,fimMetric,[],nTotalCells);
            FIMOptNextExpt = ModelGuess.totalFim(fimResults,nextExperiment+nTotalCells);

        case 'random'
            nextExperiment = randomCell(iExpt,:);
            FIMOptNextExpt = ModelGuess.totalFim(fimResults,nextExperiment+nTotalCells);

        case 'intuition'
            nextExperiment = intuitiveDesign(iExpt,:);
            FIMOptNextExpt = ModelGuess.totalFim(fimResults,nextExperiment+nTotalCells);

        case 'uniform'
            nextExperiment = uniformCell(iExpt,:);
            FIMOptNextExpt = ModelGuess.totalFim(fimResults,nextExperiment+nTotalCells);        
    end

    % Compute and Save Covariance from MH from current stage.
    parametersFound{iExpt} = newPars;
    FIMcurrentExptSaved{iExpt} = FIMCurrentExpt;
    FIMcurrentExptTrueSaved{iExpt} = FIMCurrentExpt_True;
    % FIMOptNextExptSaved{iExpt} = FIMOptNextExpt;
    covLogMH{iExpt} = cov(MHResults.mhSamples);
    covMH{iExpt} = cov(exp(MHResults.mhSamples));
    MHResultsSaved{iExpt} = MHResults;

    % Save Predicted FIM for next stage.
    for jFIM=1:nFIMsamples
        covFIM_Prediction{iExpt}{jFIM} = inv(FIMOptNextExpt{jFIM});
    end

    exptDesigns{iExpt} = nextExperiment;

    FIMOptNextExptReduced = cell(size(FIMOptNextExpt));
    for i = 1:nFIMsamples
        FIMOptNextExptReduced{i} = FIMOptNextExpt{i}(ModelGuess.fittingOptions.modelVarsToFit,...
            ModelGuess.fittingOptions.modelVarsToFit); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    if showPlots
        f = figure;
        f.Name = ['Current MH Results and Next FIM Prediction (Round ',num2str(iExpt),')'];
        ModelGuess.plotMHResults(MHResults,FIMOptNextExptReduced,fimScale,mhPlotScale,f)

        f = figure;
        f.Name = ['Current MH Results and Perfect FIM Prediction (Round ',num2str(iExpt),')'];
        ModelGuess.plotMHResults(MHResults,FIMCurrentExpt,fimScale,mhPlotScale,f)

        figure;
        title(['Number of cells Measured at each Time Point (Round ',num2str(iExpt),')']);
        bar(ModelGuess.tSpan,nTotalCells,'BarWidth', 0.4)
        ylabel('Number of Cells Measured');
        xlabel('time [min]');
        ylim([0 300]);
    end
%     if iExpt>1
%         f = figure;
%         f.Name = ['Current MH Results and Previous FIM Prediction (Round ',num2str(iExpt),')'];
%         ModelGuess.plotMHResults(MHResults,[FIMOptNextExptSaved{iExpt-1}],fimScale,mhPlotScale,f)
%     end
    saveExpName =[saveFileName,'_', data,'_',sampleType,'_',int2str(fitParameters(end))];
    save(saveExpName,'parametersFound','FIMcurrentExptSaved','FIMcurrentExptTrueSaved','covMH',...
        'covLogMH','exptDesigns','MHResultsSaved')

end


