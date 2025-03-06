classdef ExperimentDesigner

    properties
        DataFilename (1,1) string
        DataType (1,1) ExperimentalDataType
        DesignStrategy (1,1) ExperimentDesignStrategy = FIMOptimized
        FitRounds (1,1) double {mustBePositive, mustBeInteger} = 3
        Model (1,1) DiscoverableModel
        MaximumFitIterations (1,1) double {mustBePositive, mustBeInteger} = 200
        MaximumObservations (1,1) double {mustBePositive, mustBeInteger} = Inf
        MaximumRounds (1,1) double {mustBePositive, mustBeInteger} = 5
        NumberOfFIMSamples (1,1) double {mustBePositive, mustBeInteger} = 10
        ObservationQuantum (1,1) double {mustBePositive, mustBeInteger} = 1
        ObservationsCompleted (1,1) double {mustBeNonnegative, mustBeInteger}
        ObservationsPerExperiment (1,1) double {mustBePositive, mustBeInteger} = 60
        RNGSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 0
        RoundsCompleted (1,1) double {mustBeNonnegative, mustBeInteger}
        SaveExpName (1,1) string = ""
        SaveFilename (1,1) string = ""
        ShowPlots (1,1) logical = false
        TestingMode (1,1) logical = false
    end
    
    methods
        function d = ExperimentDesigner()
            
        end

        function performNextRound(designer)
            arguments
                designer (1,1) ExperimentDesigner
            end
            curRound = designer.RoundsCompleted + 1;
            clc
            disp(['Round: ',num2str(curRound)])
            nTotalCells = nTotalCells + nextExperiment;
 
    switch data
        case 'real'
            % Sample from real data
            [~,dataFile,allDataSoFar,maxAvailable] = sampleGRExperiment('ExampleData/Data.csv',...
                ModelTrue{1}.tSpan,[1,10,100],nTotalCells,curRound,datFileName);
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
            nextExperiment = randomCell(curRound,:);

        case 'intuition'
            nextExperiment = intuitiveDesign(curRound,:);

        case 'uniform'
            nextExperiment = uniformCell(curRound,:);
    end
    FIMOptNextExpt = ModelGuess{1}.totalFim(fimResults,nextExperiment+nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);

    % Compute and Save Covariance from MH from CURRENT stage.
    parametersFound{curRound} = newPars;
    covLogMH{curRound} = cov(MHResults.mhSamples);
    covMH{curRound} = cov(exp(MHResults.mhSamples));
    MHResultsSaved{curRound} = MHResults;
    
    % Save FIM predictions for NEXT stage.
    FIMcurrentExptSaved{curRound} = FIMCurrentExpt;
    FIMcurrentExptTrueSaved{curRound} = FIMCurrentExpt_True;
    FIMpredNextExpt{curRound} = FIMOptNextExpt;

    % Save Predicted COV for NEXT stage.
    for jFIM=1:nFIMsamples
        covFIM_Prediction{curRound}{jFIM} = inv(FIMOptNextExpt{jFIM});
    end

    % Reshape next experiment design
    nextExperiment = reshape(nextExperiment',[nT,nInputs])';
    exptDesigns{curRound} = nextExperiment;

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
        f.Name = ['Current MH Results and Next FIM Prediction (Round ',num2str(curRound),')'];
        ModelGuess{1}.plotMHResults(MHResults,FIMOptNextExptReduced,fimScale,mhPlotScale,f)

        if curRound>1
            f = figure;
            f.Name = ['Current MH Results and Previous FIM Prediction (Round ',num2str(curRound),')'];
            ModelGuess{1}.plotMHResults(MHResults,FIMpredNextExpt{curRound-1},fimScale,mhPlotScale,f)
        end

        f = figure;
        f.Name = ['Current MH Results and Perfect FIM Prediction (Round ',num2str(curRound),')'];
        ModelTrue{1}.plotMHResults(MHResults,FIMCurrentExpt_True,fimScale,mhPlotScale,f)

        figure;
        title(['Number of cells Measured at each Time Point (Round ',num2str(curRound),')']);
        bar(ModelGuess{1}.tSpan,nTotalCells,0.4, 'stacked')
        ylabel('Number of Cells Measured');
        xlabel('time [min]');
        ylim([0 300]);
    end

    % Save results
    save(saveExpName,'parametersFound','FIMcurrentExptSaved','FIMcurrentExptTrueSaved','covMH',...
        'covLogMH','exptDesigns','MHResultsSaved','FIMpredNextExpt')

        end

        function d = setupNamesAndPaths(d)
            arguments
                d (1,1) ExperimentDesigner
            end

            % Check if Save File Exists
            if isempty(d.SaveFilename)
                d.SaveFilename = ['results/IterativeExperimentResults_', ...
                    d.Model.ArchitectureName];
            end
            d.SaveExpName = lower([d.SaveFilename,'.mat']);
            
            % Exit if the savefile already exists.
            if exist(d.SaveExpName,'file')
                warning('Save File Already Exists -- Skipping Calculations')
                return
            else
                dirJ = strfind(d.SaveExpName,'/');
                for i=1:length(dirJ)
                    if ~exist(d.SaveExpName(1:dirJ(i)-1),"dir")
                        mkdir(d.SaveExpName(1:dirJ(i)-1));
                    end
                end
            end

            % File where data will be saved 
            d.DataFilename = d.SaveFilename;
            J = strfind(d.DataFilename,'/');
            if ~isempty(J)
                d.DataFilename= d.DataFilename(J(end)+1:end);
            end
        end % setupNamesAndPaths
    end
end

