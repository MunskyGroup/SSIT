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
else
    dirJ = strfind(saveExpName,'/');
    for i=1:length(dirJ)
        if ~exist(saveExpName(1:dirJ(i)-1),"dir")
            mkdir(saveExpName(1:dirJ(i)-1));
        end
    end

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
    case {'poisson','poissonfewer'}
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
   
    case {'burst','burstfewer'}
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

    case {'gr','grfewerinitialcells','grfewestinitialcells'}
        ModelTrue = SSIT;
        ModelTrue.species = {'cytGR';'nucGR'};
        ModelTrue.initialCondition = [20;1];
        ModelTrue.fspOptions.bounds = [0,0,30,30];
        ModelTrue.fspOptions.verbose = false;
        ModelTrue.fspOptions.fspIntegratorAbsTol = 1e-10;
        ModelTrue.sensOptions.solutionMethod = 'finiteDifference';
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
fitOptions = optimset('Display','none','MaxIter',maxFitIter);

if testing
    % When testing, use exact parameters 
    for iInput = 1:nInputs
        ModelTrue{iInput}.parameters = truePars;
        initialParameterGuess = [truePars{ModelTrue{1}.fittingOptions.modelVarsToFit,2}];
    end

    % And short inference runs.
    nSamplesMH = 500; % Number of MH Samples to run
    nFitRounds = 1;
    fitOptions.Display = 'none';
    fitOptions.MaxIter = 200;
end


ModelSolution = cell(1,nInputs);
fimTrue = cell(nInputs*nT,1);
for iInput = 1:nInputs
    ModelTrue{iInput}.fspOptions.fspTol = 1e-4;

    %% Generate Model Propensity Functions and Solve True Model
    for i=1:3
        [ModelSolution{iInput},ModelTrue{iInput}.fspOptions.bounds] = ModelTrue{iInput}.solve;
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
% Keep track of accumulated data in each experiment design
allDataSoFar = cell(1,nInputs);  
nTotalCells = zeros(nInputs,nT); 

for iExpt = 1:nExptRounds
    clc
    disp(['Round: ',num2str(iExpt)])
    nTotalCells = nTotalCells + nextExperiment;
 
    switch data
        case 'real'
            % Sample from real data
            [~,dataFile,allDataSoFar,maxAvailable] = sampleGRExperiment('ExampleData/Data.csv',...
                ModelTrue{1}.tSpan,[1,10,100],nTotalCells,iExpt,datFileName);
            dataToFit = {'cytGR','normgrcyt';'nucGR','normgrnuc'};
        case 'simulated'
            % Generate "Fake" data
            dataFile = cell(1,nInputs);
            for iInput = 1:nInputs
                % Check that tSpan has right number of time points
                if length(ModelTrue{iInput}.tSpan)~=nT
                    error('Length of tSpan does not match experiment time points.')
                end

                dataFile{iInput} = ['simData/FakeExperiment_',datFileName,'_',num2str(iInput),'.csv'];

                % Experiment settings - number of sims and 1 replica
                ModelTrue{iInput}.ssaOptions.nSimsPerExpt = max(nextExperiment(iInput,:));
                ModelTrue{iInput}.ssaOptions.Nexp = 1;
                
                if ModelTrue{iInput}.ssaOptions.nSimsPerExpt>0
                    % Sample and save the simulated data
                    ModelTrue{iInput}.sampleDataFromFSP(ModelSolution{iInput},...
                        dataFile{iInput});

                    % DownSelect fake data to specified number of cells at each time
                    X = importdata(dataFile{iInput});

                    for iT = 1:length(ModelTrue{iInput}.tSpan)
                        if nextExperiment(iInput,iT)>0
                            % Find the next nextExperiment(iInput,iT) data
                            % points at chosen time point and chosen input.
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
    stateSpaces = cell(1,nInputs);
    for iInput = 1:nInputs
        if ~isempty(allDataSoFar{iInput})&&max(sum(allDataSoFar{iInput}))>0
            hasData(iInput) = true;
            ModelGuess{iInput} = ModelGuess{iInput}.loadData(dataFile{iInput},dataToFit);
            ModelGuess{iInput}.tSpan = ModelTrue{iInput}.tSpan;
            ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
            [fspSoln,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
            % ModelGuess{iInput}.fspOptions.fspTol = inf;
            stateSpaces{iInput} = fspSoln.stateSpace;
        end
    end

    % Fit Model to data.
    for iFitIter=1:nFitRounds
        % Call function to assemble total likelihood function
        objFunc = @(x)-getObjective(x,ModelGuess,hasData,stateSpaces);
        newPars = exp(fminsearch(objFunc,log(newPars),fitOptions));
        
        % Update models with new parameters
        for iInput = 1:nInputs
            ModelGuess{iInput}.parameters(fitParameters,2) = num2cell(newPars);
            ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
            % Solve to re-set fsp bounds and state space.
            [fspSoln,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
            stateSpaces{iInput} = fspSoln.stateSpace;
            % ModelGuess{iInput}.fspOptions.fspTol = inf;
        end
    end
    
    % Show fit plots if requested.
    if showPlots
        for iInput = 1:nInputs
            if ~isempty(allDataSoFar{iInput})&&max(sum(allDataSoFar{iInput}))>0
                ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
                ModelGuess{iInput}.makeFitPlot([],1);
                % ModelGuess{iInput}.fspOptions.fspTol = inf;
            end
        end
    end

    %% MH Tuning then running.

    % Reduce the proposal distribution scale for each for each species
    % included in the fitting data.
    FIMcovScale = 1;%/(size(dataToFit,1))^2;    
    mhTune = true;
    nSamplesMHtune = 100;
    while mhTune

        %% Metropolis Hastings
        disp('Starting New MH Chain for Tuning')
        nThinMH = 2; % Thin rate for MH sampling
        nBurnMH = 1; % Number for MH burn in
        MHFitOptions.progress=false;
        MHFitOptions.saveFile = ['TMPMHChain_',datFileName,'_',num2str(ind),'.mat'];

        % Delete old file if it exists.
        delete(MHFitOptions.saveFile)

        % Call function to assemble total likelihood function
        OBJmh = @(x)getObjective(x,ModelGuess,hasData,stateSpaces);

        %% Compute FIM for use in MH Proposal Function
        if strcmp(data,'simulated') %||testing
            % Use the real FIM for the MH if we are using simulated data.
            nCellsVec = zeros(1,nInputs*nT);
            for iInput = 1:nInputs
                nCellsVec(1,((iInput-1)*nT+1:iInput*nT)) = nTotalCells(iInput,:);
            end
            FIM = ModelGuess{1}.evaluateExperiment(fimTrue,nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);

        else
            fimResults = cell(nInputs*nT,1);
            nCellsVec = zeros(1,nInputs*nT);
            for iInput = 1:nInputs
                % if ~isempty(allDataSoFar{iInput})&&sum(allDataSoFar{iInput})>0
                    ModelGuess{iInput}.tSpan = ModelTrue{iInput}.tSpan;
                    ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
                    fimResults((iInput-1)*nT+1:iInput*nT,1) =  ModelGuess{iInput}.computeFIM([],fimScale);
                    % ModelGuess{iInput}.fspOptions.fspTol = inf;
                    nCellsVec(1,((iInput-1)*nT+1:iInput*nT)) = nTotalCells(iInput,:);
                % end
            end

            % Call function to assemble full FIM from cell
            % counts and prior covariance information.
            FIM = ModelGuess{1}.evaluateExperiment(fimResults,nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);
        end
    
        FIMfree = FIM{1}(ModelGuess{1}.fittingOptions.modelVarsToFit,ModelGuess{1}.fittingOptions.modelVarsToFit);
     
        if min(eig(FIMfree))<0.1
            disp('Warning -- FIM has one or more small eigenvalues.  Reducing proposal in those directions. MH Convergence may be slow.')
            FIMfree = FIMfree + 0.1*eye(length(FIMfree));
        end
     
        % Use FIM to define MHA proposal function (FIM -> covariance of MVN)
        covFree = FIMfree^-1;
        covFree = (covFree+covFree')/2;
        proposalDistribution=@(x)mvnrnd(x,covFree * FIMcovScale);

        if testing
            MHFitOptions.progress=true;
        end

        %% Run MH Algorithm for tuning.
        [MHResults.mhSamples,MHResults.mhAcceptance,MHResults.mhValue,x0] = ...
            ssit.parest.metropolisHastingsSample(log(newPars),nSamplesMHtune,...
            'logpdf',OBJmh,'proprnd',proposalDistribution,...
            'symmetric',true,...
            'thin',nThinMH,'nchain',1,'burnin',nBurnMH,...
            'progress',MHFitOptions.progress,...
            'saveFileName',MHFitOptions.saveFile);

        delete(MHFitOptions.saveFile)
        
        mhTune = false;
        if max(abs(newPars-exp(x0))./newPars)>0  % If significant change in parameters   
            % Convert best parameters to linear space.
            newPars = exp(x0);
  
            % Update models and statespaces
            for iInput = 1:nInputs
                ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
                ModelGuess{iInput}.parameters(fitParameters,2) = num2cell(newPars);
                [fspSoln,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
                stateSpaces{iInput} = fspSoln.stateSpace;
                % ModelGuess{iInput}.fspOptions.fspTol = inf;
            end

            % Run another fminsearch starting at new parameter set.
            % Call function to assemble total likelihood function
            objFunc = @(x)-getObjective(x,ModelGuess,hasData,stateSpaces);
            newPars = exp(fminsearch(objFunc,log(newPars),fitOptions));

            % Update model parameters and statespace
            for iInput = 1:nInputs
                ModelGuess{iInput}.parameters(fitParameters,2) = num2cell(newPars);
                ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
                [fspSoln,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
                stateSpaces{iInput} = fspSoln.stateSpace;
                % ModelGuess{iInput}.fspOptions.fspTol = inf;
            end

            % Run another tuning round with new starting point.
            mhTune = true;

        elseif MHResults.mhAcceptance<0.15
            FIMcovScale = FIMcovScale/2;

            % Run another tuning round with new starting point.
            mhTune = true;

        end
    end

    disp('Starting New MH Chain')
    % Run new MH with tuned proposal.    
    [MHResults.mhSamples,MHResults.mhAcceptance,MHResults.mhValue,x0] = ...
        ssit.parest.metropolisHastingsSample(log(newPars),nSamplesMH,...
        'logpdf',OBJmh,'proprnd',proposalDistribution,...
        'symmetric',true,...
        'thin',nThinMH,'nchain',1,'burnin',nBurnMH,...
        'progress',MHFitOptions.progress,...
        'saveFileName',MHFitOptions.saveFile);

    % Update parameters if necessary.
    if max(abs(newPars-exp(x0))./newPars)>0  % If significant change in parameters
        % Convert best parameters to linear space.
        newPars = exp(x0);

        % Update models and statespaces
        for iInput = 1:nInputs
            ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
            ModelGuess{iInput}.parameters(fitParameters,2) = num2cell(newPars);
            [fspSoln,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
            stateSpaces{iInput} = fspSoln.stateSpace;
            % ModelGuess{iInput}.fspOptions.fspTol = inf;
        end
    end

    % Compute FIM for subsampling of MH results.
    J = floor(linspace(nSamplesMH/2,nSamplesMH,nFIMsamples));
    MHSamplesForFIM = exp(MHResults.mhSamples(J,:));
    
    fimResults = cell(nInputs*nT,nFIMsamples);
    for iInput = 1:nInputs
        ModelGuess{iInput}.tSpan = ModelTrue{iInput}.tSpan;
        ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
        fimResults((iInput-1)*nT+1:iInput*nT,:) =  ModelGuess{iInput}.computeFIM([],fimScale,MHSamplesForFIM);
        % ModelGuess{iInput}.fspOptions.fspTol = inf;
    end
 
    % FIM current experiment
    FIMCurrentExpt = ModelGuess{1}.totalFim(fimResults,nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);

    % True FIM for current experiment
    FIMCurrentExpt_True = ModelTrue{1}.totalFim(fimTrue,nCellsVec,ModelTrue{1}.fittingOptions.logPriorCovariance);
    
    switch lower(sampleType)
        case 'fimopt'
            % Find optimal NEXT experiment design given parameter sets
            nextExperiment = ModelGuess{1}.optimizeCellCounts(fimResults,numCellsPerExperiment,fimMetric,...
                [],nCellsVec,maxAvailable,'mean',...
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
    
    % Save FIM predictions for NEXT stage.
    FIMcurrentExptSaved{iExpt} = FIMCurrentExpt;
    FIMcurrentExptTrueSaved{iExpt} = FIMCurrentExpt_True;
    FIMpredNextExpt{iExpt} = FIMOptNextExpt;

    % Save Predicted COV for NEXT stage.
    for jFIM=1:nFIMsamples
        covFIM_Prediction{iExpt}{jFIM} = inv(FIMOptNextExpt{jFIM});
    end

    % Reshape next experiment design
    nextExperiment = reshape(nextExperiment',[nT,nInputs])';
    exptDesigns{iExpt} = nextExperiment;

    % make figures if requested.
    if showPlots
        % Plotting options
        mhPlotScale = 'log10';  % Show MH and FIM plots in log10 scale.

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


function objFun = getObjective(x,Mods,hasData,stateSpace)
% This function takes a group of models and common parameter set and then
% calculates the log-likelihood function.
arguments
    x
    Mods
    hasData
    stateSpace = [];
    % prior = @(x)0
end
if isempty(stateSpace)
    stateSpace = cell(1,length(Mods));
end
objFun = 0;
includePrior = true;
for i = 1:length(Mods)
    % Check if there is data associated with that model.
    if hasData(i)
        % We only want to include the prior in the first likelihood
        % calculation.  Otherwise it will be amplified again in every data
        % set.
        if includePrior
            objFun = objFun + Mods{i}.computeLikelihood(exp(x),stateSpace{i});
            includePrior = false;
        else
            TMP = Mods{i};
            TMP.fittingOptions.logPrior = [];
            objFun = objFun + TMP.computeLikelihood(exp(x),stateSpace{i});            
        end
    end
end
% objFun = objFun + prior(exp(x));

