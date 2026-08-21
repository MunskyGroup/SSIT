function [saveExpName] = iterativeExperimentRunnerMultiConditions(example,data,sampleType,nExptRounds,...
    rngSeed,ind,incrementAdd,numCellsPerExperiment,initialExperiment,...
    nFIMsamples,truePars,saveFileName,initialParameterGuess,...
    inputLibrary,...
    testing)
arguments
    example = 'poisson'
    data = 'simulated'
    sampleType = 'fimopt'
    nExptRounds = 5
    rngSeed = []
    ind = randi(10000)
    incrementAdd = 1
    numCellsPerExperiment = 60
    initialExperiment = []
    nFIMsamples = 10;
    truePars = [];
    saveFileName = [];
    initialParameterGuess = [];
    inputLibrary = {};
    testing=false
end

% Check if Save File Exists
if isempty(saveFileName)
    saveFileName = ['results/IterativeExperimentResults_',example];
end
saveExpName = lower([saveFileName,'.mat']);

% Exit if the savefile already exists.
if exist(saveExpName,'file')
    warning('Save File Already Exists -- Skipping Calculations')
    return
end

addpath(genpath('../../src'));
if ~isempty(rngSeed)
    rng(rngSeed)
end

% Default MLE Fitting Options
maxFitIter = 200;
nFitRounds = 3;

showPlots = false;

% Maximum number of cells is infinite unless specified otherwise.
maxAvailable = []; 


% File to save data to.
datFileName = saveFileName;
J = strfind(datFileName,'/');
if ~isempty(J)
    datFileName=datFileName(J(end)+1:end);
end


%% Define Model
switch lower(example)
    case 'poisson'
        ModelTrue = SSIT;
        ModelTrue.species = {'rna'};
        ModelTrue.initialCondition = 0;
        ModelTrue.propensityFunctions = {'kr*IDex/(kD+IDex)';'gr*rna'};
        ModelTrue.stoichiometry = [1,-1];
        if isempty(truePars)
            truePars = {'kr',10;'gr',0.3;'kD',5};
        end
        ModelTrue.parameters = truePars;
        dataToFit = {'rna','exp1_s1'};
        fitParameters = [1:3];
        nT = 21;
        ModelTrue.tSpan = linspace(0,10,nT);
        fimMetric = 'DetCovariance'; %'Determinant'';
        nSamplesMH = 1000; % Number of MH Samples to run

        %% Prior
        muLog10Prior = [0,0,0];
        sigLog10Prior = [1 1 1];
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
        ModelTrue.fittingOptions.logPriorCovariance = diag(sigLog10Prior.^2*log(10^2));

        if isempty(inputLibrary)
            ModelTrue.inputExpressions = {'IDex','5'};
            ModelTrue = ModelTrue.formPropensitiesGeneral([datFileName,'_S',num2str(ind)],true);
            ModelTrue = {ModelTrue};
        else
            TMP = cell(1,length(inputLibrary));
            for iInput = 1:length(inputLibrary)
                TMP{iInput} = ModelTrue;
                TMP{iInput}.inputExpressions = inputLibrary{iInput};
                TMP{iInput} = TMP{iInput}.formPropensitiesGeneral([datFileName,'_s',num2str(iInput)],true);
            end
            ModelTrue = TMP;
        end
   
    case 'burst'
        ModelTrue = SSIT;
        ModelTrue.species = {'on';'off';'rna'};
        ModelTrue.initialCondition = [0;1;0];
        ModelTrue.propensityFunctions = {...
            'kon*((1-2*atan(alph)/pi) + 2*atan(alph)/pi*((1e-6+IDex)/(M+((1e-6+IDex)))))*off';...
            'koff*(2*atan(alph)/pi + (1-2*atan(alph)/pi) / (((1e-6+IDex)/(M+((1e-6+IDex))))))*on';...
            'kr*on';'gr*rna'};
        ModelTrue.stoichiometry = [1,-1,0,0;...
            -1,1,0,0;...
            0,0,1,-1];

        ModelTrue.fspOptions.bounds = [0;0;0;1;1;75];

        if isempty(truePars)
            truePars = ({'kon',0.1;'koff',0.2;'kr',10;'gr',0.3;'M',4;'alph',1e-4});
        end
        ModelTrue.parameters = truePars;
        dataToFit = {'rna','exp1_s3'};
        ModelTrue.pdoOptions.unobservedSpecies = {'on','off'};

        fitParameters = [1:6];
        nT = 31;
        ModelTrue.tSpan = linspace(0,30,nT);
        fimMetric = 'GR1:5'; %'Determinant'';

        %% Prior
        muLog10Prior = [0,0,0,0,0,0];
        sigLog10Prior = [2 2 2 2 2 4];
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
        ModelTrue.fittingOptions.logPriorCovariance = diag(sigLog10Prior.^2*log(10^2));

        if isempty(inputLibrary)
            ModelTrue.inputExpressions = {'IDex','1'};
            ModelTrue = ModelTrue.formPropensitiesGeneral([datFileName,'_S',num2str(ind)],true);
            ModelTrue = {ModelTrue};
        else
            TMP = cell(1,length(inputLibrary));
            for iInput = 1:length(inputLibrary)
                TMP{iInput} = ModelTrue;
                TMP{iInput}.inputExpressions = inputLibrary{iInput};
                TMP{iInput} = TMP{iInput}.formPropensitiesGeneral([datFileName,'_s',num2str(iInput)],true);
            end
            ModelTrue = TMP;
        end
        % MLE Fitting Options
        nSamplesMH = 5000; % Number of MH Samples to run
        maxFitIter = 2000;
        nFitRounds = 5;

    case 'gr'
        ModelTrue = SSIT;
        ModelTrue.species = {'cytGR';'nucGR'};
        ModelTrue.initialCondition = [20;1];
        ModelTrue.fspOptions.bounds = [0,0,30,30];
        ModelTrue.fspOptions.verbose = false;
        ModelTrue.fspOptions.fspIntegratorAbsTol = 1e-10;
        ModelTrue.propensityFunctions = {'kcn0*(1 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';...
            'knc*nucGR'; 'kg1';'gnuc*nucGR'};
        ModelTrue.stoichiometry = [-1,1,1,0;...
            1,-1,0,-1];
        ModelTrue.customConstraintFuns = {'x1+x2'};

        % Provide parameter guess -- not used, but needed to create model.
        if isempty(truePars)
            truePars = {'kcn0',0.005;'kcn1',0.08;'knc',0.014;...
                'kg1',0.012;'gnuc',0.005;'MDex',10.44};
            ModelTrue.parameters = truePars;
        end
        
        ModelTrue.fspOptions.initApproxSS = true;
        ModelTrue.fspOptions.usePiecewiseFSP = true;
        ModelTrue.fspOptions.constantJacobian = true;
        ModelTrue.fspOptions.constantJacobianTime = 0.1;

        dataToFit = {'cytGR','exp1_s1';'nucGR','exp1_s2'};
        % dataToFit = {'cytGR','exp1_s1'};
        fitParameters = [1:6];
        nT = 6;
        ModelTrue.tSpan = [0,10,30,50,75,180];
        fimMetric = 'GR[2,5,6]'; %'Determinant'';

        %% Prior
        muLog10Prior = [-2 1 -2 -2 -2 1];
        sigLog10Prior = 2*ones(1,6);
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
        ModelTrue.fittingOptions.logPriorCovariance = diag(sigLog10Prior.^2*log(10^2));

        % Create Separate Models for each Provided Input Signal
        if isempty(inputLibrary)
            ModelTrue.inputExpressions = {'IDex','100'};
            ModelTrue = ModelTrue.formPropensitiesGeneral([datFileName,'_S',num2str(ind)],true);
            ModelTrue = {ModelTrue};
        else
            TMP = cell(1,length(inputLibrary));
            for iInput = 1:length(inputLibrary)
                TMP{iInput} = ModelTrue;
                TMP{iInput}.inputExpressions = inputLibrary{iInput};
                TMP{iInput} = TMP{iInput}.formPropensitiesGeneral([datFileName,'_s',num2str(iInput)],true);
            end
            ModelTrue = TMP;
        end

        % MLE Fitting Options
        nSamplesMH = 5000; % Number of MH Samples to run
        maxFitIter = 1000;
        nFitRounds = 5;

end
% Shut down the parallel pool to do rest of the work on the CPU.
% delete(gcp('nocreate'));

nInputs = length(inputLibrary);
fitOptions = optimset('Display','iter','MaxIter',maxFitIter);

if testing
    % When testing, use exact parameters 
    for iInput = 1:nInputs
        ModelTrue{iInput}.parameters = truePars;
        initialParameterGuess = [truePars{ModelTrue{1}.fittingOptions.modelVarsToFit,2}];
    end

    % And short inference runs.
    nSamplesMH = 500; % Number of MH Samples to run
    nFitRounds = 1;
    fitOptions.Display = 'iter';
    fitOptions.MaxIter = 200;
end


ModelSolution = cell(1,nInputs);
fimTrue = cell(nInputs*nT,1);
for iInput = 1:nInputs
    ModelTrue{iInput}.fspOptions.fspTol = 1e-4;

    %% Generate Model Propensity Functions and Solve True Model
    for i=1:3
        ModelTrue{iInput} = ModelTrue{iInput}.solve(solver = "fsp");
        %[ModelSolution{iInput},ModelTrue{iInput}.fspOptions.bounds] = ModelTrue{iInput}.solve;
    end

    %% FIM options
    fimScale = 'log'; % Maximize fim for log parameters

    %% True Model FIM
    ModelTrue{iInput}.fspOptions.fspTol = 1e-8;
    fimTrue((iInput-1)*nT+1:iInput*nT) = ModelTrue{iInput}.computeFIM([],fimScale);

    %% Verify that the true model and simulated data look correct.
    if showPlots
        dataFile = ['simData/FakeExperiment_',datFileName,'_',num2str(iInput),'.csv'];
        nextExperiment = 100*ones(1,nT);
        ModelTrue{iInput}.ssaOptions.nSimsPerExpt = max(nextExperiment);
        ModelTrue{iInput}.ssaOptions.Nexp = 1;
        ModelTrue{iInput}.sampleDataFromFSP(ModelSolution{iInput},dataFile);
        ModelTMP = ModelTrue{iInput}.loadData(dataFile,dataToFit);
        ModelTMP.makeFitPlot([],1);
    end
end

%% Set Rules for Experiment Designs and Iterations
% Possible Time points for Experiments.
    
% Definition of initial experiment
if isempty(initialExperiment)
    initialExperiment = zeros(nInputs,nT);
    initialExperiment(1,[1,round(nT/3),round(2*nT/3),nT]) = round(numCellsPerExperiment/4);
end

nextExperiment = initialExperiment;

%% Experiment Design Definitions
% Random distribution of experiments
N = numCellsPerExperiment/incrementAdd; % needs to be int
randomCell = zeros(nExptRounds,nInputs*nT);
for i = 1:nExptRounds
    n = randi(nT*nInputs,1,N);
    for j = n
        randomCell(i,j) = randomCell(i,j) + incrementAdd;
    end
end

% Uniform distribution of experiments
uniformCell = floor(ones(nExptRounds,nInputs*nT)*numCellsPerExperiment/(nInputs*nT));
for i=1:nExptRounds
    if sum(uniformCell(i,:))<numCellsPerExperiment
        uniformCell(i,1:numCellsPerExperiment-sum(uniformCell(i,:)))=...
            uniformCell(i,1:numCellsPerExperiment-sum(uniformCell(i,:)))+1;
    end
end

%% Create Model for Estimate.
ModelGuess = cell(1,nInputs);
for iInput = 1:nInputs
    ModelTrue{iInput}.fittingOptions.modelVarsToFit = fitParameters;
    ModelGuess{iInput} = ModelTrue{iInput};
end