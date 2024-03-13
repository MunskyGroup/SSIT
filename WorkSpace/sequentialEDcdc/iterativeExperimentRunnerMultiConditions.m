function [saveExpName] = iterativeExperimentRunnerMultiConditions(example,data,sampleType,nExptRounds,...
    rngSeed,ind,incrementAdd,numCellsPerExperiment,initialExperiment,...
    nFIMsamples,truePars,saveFileName,initialParameterGuess,...
    inputLibrary)
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
end

% Check if Save File Exists
if isempty(saveFileName)
    saveFileName = ['IterativeExperimentResults_',example];
end
saveExpName = lower([saveFileName,'_', data,'_',sampleType,'_',int2str(ind),'.mat']);
if exist(saveExpName,'file')
    warning('Save File Already Exists -- Skipping Calculations')
    return
end

addpath(genpath('../../src'));
if ~isempty(rngSeed)
    rng(rngSeed)
end

showPlots = false;

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
        fimMetric = 'Determinant';
        nSamplesMH = 1000; % Number of MH Samples to run
        
        %% Prior
        muLog10Prior = [0,0,0];
        sigLog10Prior = [1 1 1];
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
        ModelTrue.fittingOptions.logPriorCovariance = diag(sigLog10Prior.^2*log(10^2));

        if isempty(inputLibrary)
            ModelTrue.inputExpressions = {'IDex','5'};
            ModelTrue = ModelTrue.formPropensitiesGeneral(['poisson_',num2str(ind)],true);
            ModelTrue = {ModelTrue};
        else
            TMP = cell(1,length(inputLibrary));
            for iInput = 1:length(inputLibrary)
                TMP{iInput} = ModelTrue;
                TMP{iInput}.inputExpressions = inputLibrary{iInput};
                TMP{iInput} = TMP{iInput}.formPropensitiesGeneral(['poisson_',num2str(iInput)],true);
            end
            ModelTrue = TMP;
        end


    case 'burst'
        ModelTrue = SSIT;
        ModelTrue.species = {'on';'off';'rna'};
        ModelTrue.initialCondition = [0;1;0];
        ModelTrue.propensityFunctions = {'kon*off';'koff*on';'kr*on';'gr*rna'};
        ModelTrue.stoichiometry = [1,-1,0,0;...
            -1,1,0,0;...
            0,0,1,-1];
        if isempty(truePars)
            truePars = ({'kon',0.1;'koff',0.2;'kr',10;'gr',0.3});
        end
        ModelTrue.parameters = truePars;
        ModelTrue = ModelTrue.formPropensitiesGeneral(['Burst_',num2str(ind)],true);
        dataToFit = {'rna','exp1_s3'};
        fitParameters = [1:4];
        nT = 61;
        ModelTrue.tSpan = linspace(0,30,nT);
        fimMetric = 'Determinant';
        nSamplesMH = 5000; % Number of MH Samples to run

        %% Prior
        muLog10Prior = [0,0,0,0];
        sigLog10Prior = [2 2 2 2];
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
        ModelTrue.fittingOptions.logPriorCovariance = diag(sigLog10Prior.^2*log(10^2));

    case 'uncertainburst'
        ModelTrue = SSIT;
        ModelTrue.species = {'on';'off';'rna'};
        ModelTrue.initialCondition = [0;1;0];
        ModelTrue.propensityFunctions = {'kon*((1-2*atan(alph)/pi)+2*atan(alph)/pi*Iu_on)*off';...
            'koff*(2*atan(alph)/pi+(1-2*atan(alph)/pi)*Iu_off)*on';...
            'kr*on';'gr*rna'};
        ModelTrue.stoichiometry = [1,-1,0,0;...
            -1,1,0,0;...
            0,0,1,-1];
        
        ModelTrue.inputExpressions = {'Iu_on','(t>0)*(t^2/(10^2+t^2))';...
            'Iu_off','(t>0)*(10^2/(10^2+t^2))'};
        
        if isempty(truePars)
            truePars = ({'kon',0.1;'koff',0.2;'kr',10;'gr',0.3;'alph',1e-4});
        end
        
        ModelTrue.parameters = truePars;
        ModelTrue = ModelTrue.formPropensitiesGeneral(saveFileName,true);
        dataToFit = {'rna','exp1_s3'};
        fitParameters = [1:5];
        nT = 61;
        ModelTrue.tSpan = linspace(0,30,nT);
        fimMetric = 'Determinant';
        nSamplesMH = 1000; % Number of MH Samples to run

        %% Prior
        muLog10Prior = [0,0,0,0,0];
        sigLog10Prior = [2 2 2 2 2];
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
        ModelTrue.fittingOptions.logPriorCovariance = diag(sigLog10Prior.^2*log(10^2));
    
    case 'toggle'
        ModelTrue = SSIT;
        if isempty(truePars)
            truePars = {'kb',10;'ka',80;'M',20;'g',1};
        end
        ModelTrue.parameters = truePars;
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
        
        nSamplesMH = 5000; % Number of MH Samples to run

        ModelTrue = ModelTrue.formPropensitiesGeneral(['Toggle_',num2str(ind)],true);

    case 'dusp1'
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
        
        %ModelTrue.parameters = ({'koff',0.14;'kon',1;'kr',1;'gr',0.01;...
        %                   'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
        if isempty(truePars)
            truePars = load('SGRS_model_v1.mat').SGRS_Model.parameters;
        end
        ModelTrue.parameters = truePars;     
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
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        muLog10Prior = muLog10Prior(ModelTrue.fittingOptions.modelVarsToFit);
        sigLog10Prior = sigLog10Prior(ModelTrue.fittingOptions.modelVarsToFit);
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);

        nSamplesMH = 5000; % Number of MH Samples to run

%         intuitiveDesign = [50 0 50 0 50 0 50 0 50 0 0 50;
%                              75 0 0 75 0 0 75 0 0 0 75 0;
%                              0 0 0 100 0 0 100 0 0 100 0 0;
%                              0 0 0 100 0 0 0 100 0 100 0 0;
%                              0 0 0 0 100 0 100 0 0 100 0 0];

                    intuitiveDesign = [20 0 20 0 20 0 20 0 20 0 0 20;
                                     30 0 0 30 0 0 30 0 0 0 30 0;
                                     0 0 0 40 0 0 40 0 0 40 0 0;
                                     0 0 0 40 0 0 0 40 0 40 0 0;
                                     0 0 0 0 40 0 40 0 0 40 0 0;
                                     0 0 40 0 40 0 0 40 0 0 0 0;
                                     40 0 0 0 40 0 0 0 0 0 0 0;
                                     0 0 0 40 0 0 0 40 0 40 0 0;
                                     0 0 0 0 40 0 40 0 0 40 0 0;
                                     80 0 0 0 40 0 0 0 0 0 0 0];

        ModelTrue = ModelTrue.formPropensitiesGeneral(['DUSP1_',num2str(ind)],true);

    case 'dexgr'
        ModelTrue = SSIT;
        ModelTrue.species = {'cytGR';'nucGR'};
        ModelTrue.initialCondition = [20;1];
        ModelTrue.fspOptions.bounds = [0,0,30,30];
        ModelTrue.fspOptions.verbose = true;
        ModelTrue.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';...
            'kg1';'gnuc*nucGR'};
        ModelTrue.stoichiometry = [-1,1,1,0;...
            1,-1,0,-1];

        if isempty(truePars)
            truePars = {'kcn0',0.005;'kcn1',0.08;'knc',0.014;...
                'kg1',0.012;'gnuc',0.005;'MDex',10.44};
        end

        ModelTrue.parameters = truePars;

        ModelTrue.fspOptions.initApproxSS = true;
        ModelTrue.fspOptions.usePiecewiseFSP = true;
        ModelTrue.fspOptions.constantJacobian = true;
        ModelTrue.fspOptions.constantJacobianTime = 1.1;

        dataToFit = {'cytGR','exp1_s1';'nucGR','exp1_s2'};
        fitParameters = [1:6];
        nT = 7;
        ModelTrue.tSpan = [0,10,30,50,75,150,180];
        fimMetric = 'Determinant';
        nSamplesMH = 1000; % Number of MH Samples to run

        %% Prior
        muLog10Prior = [-2 -1 -1 -2 -2 1];
        sigLog10Prior = 2*ones(1,6);
        ModelTrue.fittingOptions.modelVarsToFit = fitParameters;
        ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
        ModelTrue.fittingOptions.logPriorCovariance = diag(sigLog10Prior.^2*log(10^2));
       
        if isempty(inputLibrary)
            ModelTrue.inputExpressions = {'IDex','100'};
            ModelTrue = ModelTrue.formPropensitiesGeneral(['grdex_',num2str(ind)],true);
            ModelTrue = {ModelTrue};
        else
            TMP = cell(1,length(inputLibrary));
            for iInput = 1:length(inputLibrary)
                TMP{iInput} = ModelTrue;
                TMP{iInput}.inputExpressions = inputLibrary{iInput};
                TMP{iInput} = TMP{iInput}.formPropensitiesGeneral(['grdex_',num2str(iInput)],true);
            end
            ModelTrue = TMP;
        end

end

nInputs = length(inputLibrary);
ModelSolution = cell(1,nInputs);
fimTrue = {};
for iInput = 1:nInputs
    ModelTrue{iInput}.fspOptions.fspTol = 1e-6;

    %% Generate Model Propensity Functions and Solve True Model
    [ModelSolution{iInput},ModelTrue{iInput}.fspOptions.bounds] = ModelTrue{iInput}.solve;

    %% FIM options
    fimScale = 'log'; % Maximize fim for log parameters

    %% True Model FIM
    fimTrue = [fimTrue;ModelTrue{iInput}.computeFIM([],fimScale)];

    %% Verify that the true model and simulated data look correct.
    if showPlots
        dataFile = ['FakeExperiment_',saveFileName,'_',num2str(iInput),'.csv'];
        nextExperiment = 100*ones(1,nT);
        ModelTrue{iInput}.ssaOptions.nSimsPerExpt = max(nextExperiment);
        ModelTrue{iInput}.ssaOptions.Nexp = 1;
        ModelTrue{iInput}.sampleDataFromFSP(ModelSolution{iInput},dataFile)
        ModelTMP = ModelTrue{iInput}.loadData(dataFile,dataToFit);
        ModelTMP.makeFitPlot([],1);
    end
end

%% Set Rules for Experiment Designs and Iterations
% Possible Time points for Experiments.
    
% MLE Fitting Options
maxFitIter = 200;
nFitRounds = 3;

% Metropolis Hastings Properties
nThinMH = 2; % Thin rate for MH sampling
nBurnMH = 100; % Number for MH burn in

% Definition of initial experiment
if isempty(initialExperiment)
    initialExperiment = zeros(nInputs,nT);
    initialExperiment(1,[1,round(nT/3),round(2*nT/3),nT]) = round(numCellsPerExperiment/4);
end

nextExperiment = initialExperiment;
% Plotting options
mhPlotScale = 'log10';  % Show MH and FIM plots in log10 scale.

%% Experiment Design Definitions
N = numCellsPerExperiment/incrementAdd; % needs to be int
randomCell = zeros(nExptRounds,nInputs*nT);
for i = 1:nExptRounds
    n = randi(nT*nInputs,1,N);
    for j = n
        randomCell(i,j) = randomCell(i,j) + incrementAdd;
    end
end

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
% Randomize the initial parameter set.
newPars = initialParameterGuess;

%% Initialize cells for saving results
covMH = cell(1,nExptRounds);
covLogMH = cell(1,nExptRounds);
covLogFIM_Prediction = cell(1,nExptRounds);
covFIM_Prediction = cell(1,nExptRounds);
parametersFound = cell(1,nExptRounds);
exptDesigns = cell(1,nExptRounds);
FIMpredNextExpt = cell(1,nExptRounds);
FIMcurrentExptSaved = cell(1,nExptRounds);
FIMcurrentExptTrueSaved = cell(1,nExptRounds);
MHResultsSaved = cell(1,nExptRounds);

%% Loop over subsequent experiments
allDataSoFar = cell(1,nInputs);
nTotalCells = zeros(nInputs,nT);

for iExpt = 1:nExptRounds
    clc
    disp(['Round: ',num2str(iExpt)])
    nTotalCells = nTotalCells + nextExperiment;
 
    switch data
        case 'real'
            % Sample from real data
            [simData,dataFile,allDataSoFar] = sampleGRExperiment('ExampleData/Gated_dataframe_Ron_030624_NormalizedGR_bins.csv',...
                ModelTrue{1}.tSpan,[1,10,100],nTotalCells,iExpt);
            dataToFit = {'cytGR','normgrcyt';'nucGR','normgrnuc'};
        case 'simulated'
            % Generate "Fake" data
            dataFile = cell(1,nInputs);
            for iInput = 1:nInputs
                dataFile{iInput} = ['FakeExperiment_',saveFileName,'_',num2str(iInput),'.csv'];

                ModelTrue{iInput}.ssaOptions.nSimsPerExpt = max(nextExperiment(iInput,:));
                ModelTrue{iInput}.ssaOptions.Nexp = 1;
                if ModelTrue{iInput}.ssaOptions.nSimsPerExpt>0
                    ModelTrue{iInput}.sampleDataFromFSP(ModelSolution{iInput},dataFile{iInput})

                    % DownSelect fake data to specified number of cells at each time
                    X = importdata(dataFile{iInput});
                    for iT = 1:length(ModelTrue{iInput}.tSpan)
                        if nextExperiment(iInput,iT)>0
                            J = find(X.data(:,1)==ModelTrue{iInput}.tSpan(iT),nextExperiment(iInput,iT));

                            % Add new data to previous data set
                            allDataSoFar{iInput} = [allDataSoFar{iInput};X.data(J,:)];
                        end
                    end
                end
                % Save generated data to file.
                A = table;
                if ~isempty(allDataSoFar{iInput})
                    for j=1:length(X.colheaders)
                        A.(X.colheaders{j}) = allDataSoFar{iInput}(:,j);
                        writetable(A,dataFile{iInput})
                    end
                end
            end

    end

    % Add Data to Model
    hasData = zeros(1,nInputs,'logical');
    for iInput = 1:nInputs
        if ~isempty(allDataSoFar{iInput})
            hasData(iInput) = true;
            ModelGuess{iInput} = ModelGuess{iInput}.loadData(dataFile{iInput},dataToFit);
            ModelGuess{iInput}.tSpan = ModelTrue{iInput}.tSpan;
            ModelGuess{iInput}.fspOptions.fspTol = inf;
        end
    end

    % Call function to assemble total likelihood function
    objFunc = @(x)-getObjective(x,ModelGuess,hasData);

    % Fit Model to data.
    fitOptions = optimset('Display','none','MaxIter',maxFitIter);
    for iFitIter=1:nFitRounds
        newPars = exp(fminsearch(objFunc,log(newPars),fitOptions));
        for i = 1:nInputs
            ModelGuess{iInput}.parameters(fitParameters,2) = num2cell(newPars);
            ModelGuess{iInput}.fspOptions.fspTol = 1e-6;
            [~,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
        end
    end
    
    % Show fit plotsif requested.
    if showPlots
        for iInput = 1:nInputs
            if ~isempty(allDataSoFar{iInput})
                ModelGuess{iInput}.makeFitPlot([],1);
            end
        end
    end

   % Run Metropolis Hastings
    disp('Starting MH Chain')
    MHFitOptions.progress=false;
    MHFitOptions.saveFile = ['TMPMHChain_',saveFileName,'_',num2str(ind),'.mat'];
    delete(MHFitOptions.saveFile) 

    % Call function to assemble total likelihood function
    for iInput = 1:nInputs
        ModelGuess{iInput}.fspOptions.fspTol = inf;
    end
    OBJmh = @(x)getObjective(x,ModelGuess,hasData);

    %% Compute FIM for use in MH Proposal Function

    fimResults = {};
    nCellsVec = [];
    for iInput = 1:nInputs
        % if ~isempty(allDataSoFar{iInput})
            fimResults =  [fimResults;ModelGuess{iInput}.computeFIM([],'log')];
            nCellsVec = [nCellsVec,nTotalCells(iInput,:)];
        % end
    end

    % % Call function to assemble full FIM from cell
    % % counts and prior covariance information.
    FIM = ModelGuess{1}.evaluateExperiment(fimResults,nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);
    % 
    FIMfree = FIM{1}(ModelGuess{1}.fittingOptions.modelVarsToFit,ModelGuess{1}.fittingOptions.modelVarsToFit);
    % 
    if min(eig(FIMfree))<1
        disp('Warning -- FIM has one or more small eigenvalues.  Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
        FIMfree = FIMfree + 1*eye(length(FIMfree));
    end
    % 
    covFree = FIMfree^-1;
    covFree = (covFree+covFree')/2;
    proposalDistribution=@(x)mvnrnd(x,covFree);

    %% Run MH Algorithm
    [MHResults.mhSamples,MHResults.mhAcceptance,MHResults.mhValue,x0] = ...
        ssit.parest.metropolisHastingsSample(log(newPars),nSamplesMH,...
        'logpdf',OBJmh,'proprnd',proposalDistribution,...
        'symmetric',true,...
        'thin',nThinMH,'nchain',1,'burnin',nBurnMH,...
        'progress',MHFitOptions.progress,...
        'saveFileName',MHFitOptions.saveFile);

    newPars = exp(x0);
  
    for iInputs = 1:nInputs
        ModelGuess{iInput}.fspOptions.fspTol = 1e-6;
        ModelGuess{iInput}.parameters(ModelGuess{iInput}.fittingOptions.modelVarsToFit,2) = num2cell(newPars);
        [~,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
    end
    delete(MHFitOptions.saveFile)

    % Compute FIM for subsampling of MH results.
    J = floor(linspace(nSamplesMH/2,nSamplesMH,nFIMsamples));
    MHSamplesForFIM = exp(MHResults.mhSamples(J,:));
    
    fimResults = {};
    for iInput = 1:nInputs
        fimResults =  [fimResults;ModelGuess{iInput}.computeFIM([],fimScale,MHSamplesForFIM)];
    end
 
    % FIM current experiment
    FIMCurrentExpt = ModelGuess{1}.totalFim(fimResults,nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);

    % True FIM for chosen experiment
    FIMCurrentExpt_True = ModelTrue{1}.totalFim(fimTrue,nCellsVec,ModelTrue{1}.fittingOptions.logPriorCovariance);
    
    switch lower(sampleType)
        case 'fimopt'
            % Find optimal NEXT experiment design given parameter sets
            nextExperiment = ModelGuess{1}.optimizeCellCounts(fimResults,numCellsPerExperiment,fimMetric,...
                [],nCellsVec,[],'mean',...
                ModelGuess{1}.fittingOptions.logPriorCovariance,...
                incrementAdd);
            
        case 'random'
            nextExperiment = randomCell(iExpt,:);

        case 'intuition'
            nextExperiment = intuitiveDesign(iExpt,:);

        case 'uniform'
            nextExperiment = uniformCell(iExpt,:);
    end
    FIMOptNextExpt = ModelGuess{1}.totalFim(fimResults,nextExperiment+nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);

    % Compute and Save Covariance from MH from CURRENT stage.
    parametersFound{iExpt} = newPars;
    covLogMH{iExpt} = cov(MHResults.mhSamples);
    covMH{iExpt} = cov(exp(MHResults.mhSamples));
    MHResultsSaved{iExpt} = MHResults;
    
    % Save FIM predictions for 
    FIMcurrentExptSaved{iExpt} = FIMCurrentExpt;
    FIMcurrentExptTrueSaved{iExpt} = FIMCurrentExpt_True;
    FIMpredNextExpt{iExpt} = FIMOptNextExpt;

    % Save Predicted FIM for next stage.
    for jFIM=1:nFIMsamples
        covFIM_Prediction{iExpt}{jFIM} = inv(FIMOptNextExpt{jFIM});
    end

    nextExperiment = reshape(nextExperiment',[nT,nInputs])';
    exptDesigns{iExpt} = nextExperiment;

    % make figures if requested.
    if showPlots
        FIMOptNextExptReduced = cell(size(FIMOptNextExpt));
        for i = 1:nFIMsamples
            FIMOptNextExptReduced{i} = FIMOptNextExpt{i}(ModelGuess{1}.fittingOptions.modelVarsToFit,...
                ModelGuess{1}.fittingOptions.modelVarsToFit); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        f = figure;
        f.Name = ['Current MH Results and Next FIM Prediction (Round ',num2str(iExpt),')'];
        ModelGuess{1}.plotMHResults(MHResults,FIMOptNextExptReduced,fimScale,mhPlotScale,f)

        if iExpt>1
            f = figure;
            f.Name = ['Current MH Results and Previous FIM Prediction (Round ',num2str(iExpt),')'];
            ModelGuess{1}.plotMHResults(MHResults,FIMpredNextExpt{iExpt-1},fimScale,mhPlotScale,f)
        end

        f = figure;
        f.Name = ['Current MH Results and Perfect FIM Prediction (Round ',num2str(iExpt),')'];
        ModelTrue{1}.plotMHResults(MHResults,FIMCurrentExpt_True,fimScale,mhPlotScale,f)

        figure;
        title(['Number of cells Measured at each Time Point (Round ',num2str(iExpt),')']);
        bar(ModelGuess{1}.tSpan,nTotalCells,0.4, 'stacked')
        ylabel('Number of Cells Measured');
        xlabel('time [min]');
        ylim([0 300]);
    end

    % Save results
    save(saveExpName,'parametersFound','FIMcurrentExptSaved','FIMcurrentExptTrueSaved','covMH',...
        'covLogMH','exptDesigns','MHResultsSaved','FIMpredNextExpt')

end


function objFun = getObjective(x,Mods,hasData,stateSpace,prior)
arguments
    x
    Mods
    hasData
    stateSpace = []
    prior = @(x)0
end
objFun = 0;
for i = 1:length(Mods)
    if hasData(i)
        objFun = objFun + Mods{i}.computeLikelihood(exp(x),stateSpace);
    end
end
objFun = objFun + prior(exp(x));

