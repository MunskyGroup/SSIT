% Gene Expression System FIM-Based SED Controller Class Definition
% This calculates the "uncertainty volume" (the determinant of the
% covariance matrix) for the posterior distribution after each experiment
% and selects the next experiment such that the expected subsequent
% uncertainty volume is minimized.
classdef GeneExprSysFIMBasedSEDController < GeneExprSysSEDController
    properties (SetAccess = protected)
        ExistingAndSimulatedDataSoFar (1,:) cell
        cellDistributions (1,1) GeneExprSysCellDistribution
        covMH (1,:) cell
        covLogMH (1,:) cell
        covLogFIM_Prediction (1,:) cell
        covFIM_Prediction (1,:) cell
        exptDesigns (1,:) cell
        FIMpredNextExpt (1,:) cell
        FIMcurrentExptSaved (1,:) cell
        FIMcurrentExptTrueSaved (1,:) cell
        FIMsamples (1,1) integer {mustBePositive} = 10
        MHResultsSaved (1,:) cell
        nInputs (1,1) integer {mustBePositive}
        nT (1,1) integer {mustBePositive}
        parametersFound (1,:) cell
        showPlots (1,1) boolean = false
    end

    methods
        % Need constructor to allow for extra properties
        function c = GeneExprSysFIMBasedSEDController(...
                model, rounds, FIMsamples, showPlots)
            superargs = {model, rounds};
            % Call superclass constructor
            c@GeneExprSysSEDController(superargs{:});

            if (nargin > 2)
                c.FIMsamples = FIMsamples;
            end
            if (nargin > 3)
                c.showPlots = showPlots;
            end

            % Initialize cell arrays for saving results
            c.covMH = cell(1,rounds);
            c.covLogMH = cell(1,rounds);
            c.covLogFIM_Prediction = cell(1,rounds);
            c.covFIM_Prediction = cell(1,rounds);
            c.parametersFound = cell(1,rounds);
            c.exptDesigns = cell(1,rounds);
            c.FIMpredNextExpt = cell(1,rounds);
            c.FIMcurrentExptSaved = cell(1,rounds);
            c.FIMcurrentExptTrueSaved = cell(1,rounds);
            c.MHResultsSaved = cell(1,rounds);

            % Initialize cell distribution matrix for possible experiments
            c.nT = length(c.EstimatedModel.tSpan);
            c.nInputs = length(c.InputLibrary);
            c.cellDistributions = GeneExprSysCellDistribution(...
                c.nInputs, c.nT);
            c.cellDistributions.Matrix(...
                1, [1, round(nT/3), round(2*nT/3), nT]) = ...
                round(c.NumCellsPerExperiment/4);
        end
    end

    methods(Access = ?GeneExprSysController)
        function experiment = SelectNextExperiment(controller)
            %% 1. Create data structures for all possible inputs:
            ModelGuess = cell(1, controller.nInputs);
            dataFile = cell(1, controller.nInputs); % Holds simulated data
            hasData = zeros(1, controller.nInputs, 'logical');
            stateSpaces = cell(1, controller.nInputs);
            nextRound = controller.CompletedRounds + 1;         

            %% 2. Evaluate each possible input:
            for iInput = 1:controller.nInputs
                % 2.a. Set up the guessed model:

                % Clone the model
                ModelGuess{iInput} = controller.ExpectedModel;
                % Adjust the span
                ModelGuess{iInput}.tSpan = newTimeSpan;
                % Adjust the input expression
                ModelGuess{iInput}.inputExpressions = controller.inputLibrary{iInput};

                % 2.b. Generate simulated data:
                
                % Destination file for data:
                dataFile{iInput} = [...
                    'simData/Experiment_', saveFileName, ...
                    'Round_', num2str(nextRound), ...
                    'Input_', num2str(iInput), '.csv'];

                % Experiment settings - number of sims and 1 replica
                ModelGuess{iInput}.ssaOptions.nSimsPerExpt = ...
                    max(controller.cellDistributions.Matrix(iInput, :));
                ModelGuess{iInput}.ssaOptions.Nexp = 1;
                
                if ModelGuess{iInput}.ssaOptions.nSimsPerExpt > 0
                    % Sample and save the simulated data
                    ModelGuess{iInput}.sampleDataFromFSP(...
                        ModelGuess{iInput}, dataFile{iInput})

                    % Downselect simulated data to specified number of
                    % cells at the next time
                    X = importdata(dataFile{iInput});

                    numCellsForNextExperiment = ...
                        controller.cellDistributions.Matrix(iInput, nextRound);
                    if numCellsForNextExperiment > 0
                        % Find the data points at the next time and chosen
                        % input:
                        J = find(X.data(:,1) == ...
                            ModelGuess{iInput}.tSpan(1),...
                            numCellsForNextExperiment);

                        % Add new data to previous data set:
                        controller.ExistingAndSimulatedDataSoFar{iInput} = ...
                            [controller.AllDataSoFar{iInput}; X.data(J,:)];
                    end
                end % if ModelGuess{iInput}.ssaOptions.nSimsPerExpt > 0

                % 2.c. Save generated data to file:
                A = table;
                if ~isempty(controller.ExistingAndSimulatedDataSoFar{iInput})
                    for j=1:length(X.colheaders)
                        A.(X.colheaders{j}) = ...
                            controller.ExistingAndSimulatedDataSoFar{iInput}(:, j);
                        writetable(A, dataFile{iInput})
                    end
                end

                % 2.d. Add Data to Model:

                if ~isempty(controller.ExistingAndSimulatedDataSoFar{iInput}) ...
                        &&all(sum(controller.ExistingAndSimulatedDataSoFar{iInput})>0)
                    hasData(iInput) = true;
                    ModelGuess{iInput} = ModelGuess{iInput}.loadData(...
                        dataFile{iInput}, dataToFit);
                    ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
                    [fspSoln, ModelGuess{iInput}.fspOptions.bounds] = ...
                        ModelGuess{iInput}.solve;
                    stateSpaces{iInput} = fspSoln.stateSpace;
                end

                % 2.e. Fit Model to data:

                for iFitIter=1:nFitRounds
                    % Call function to assemble total likelihood function
                    objFunc = @(x)-getObjective(x, ModelGuess, hasData, stateSpaces);
                    newPars = exp(fminsearch(objFunc, log(newPars), fitOptions));
                    
                    % Update models with new parameters
                    
                    ModelGuess{iInput}.parameters(fitParameters, 2) = ...
                        num2cell(newPars);
                    ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
                    % Solve to re-set fsp bounds and state space.
                    [fspSoln, ModelGuess{iInput}.fspOptions.bounds] = ...
                        ModelGuess{iInput}.solve;
                    stateSpaces{iInput} = fspSoln.stateSpace;                      
                end

                % 2.f. Show fit plots if requested:
                if controller.showPlots                   
                    if ~isempty(controller.ExistingAndSimulatedDataSoFar{iInput}) ...
                            &&sum(controller.ExistingAndSimulatedDataSoFar{iInput})>0
                        ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
                        ModelGuess{iInput}.makeFitPlot([], 1);
                    end                   
                end

                % 2.g. Tune and run the Metropolis-Hastings algorithm:

                controller.TuneAndRunMH(iInput); 
            end % for iInput ...
        end

        function updateModel(controller, data)
        end
    end

    methods(Access = private)
        function TuneAndRunMH(controller, iInput, ...
                ModelGuess, dataFile, hasData, stateSpaces)
            arguments
                controller (1,1) GeneExprSysSEDController
                iInput (1,1) integer {mustBePositive}
                ModelGuess (1,:) cell
                dataFile (1,:) cell
                hasData (1,:) logical
                stateSpaces (1,:) cell
            end
            
            % Reduce the proposal distribution scale for each species
            % included in the fitting data.
            FIMcovScale = 1;
            mhTune = true;
            nSamplesMHtune = 100;
            while mhTune
        
                %% Metropolis Hastings
                disp('Starting New MH Chain for Tuning')
                nThinMH = 2; % Thin rate for MH sampling
                nBurnMH = 1; % Number for MH burn in
                MHFitOptions.progress=false;
                MHFitOptions.saveFile = ...
                    ['TMPMHChain_', saveFileName, '_', num2str(iInput), '.mat'];
        
                % Delete old file if it exists.
                delete(MHFitOptions.saveFile)
        
                % Call function to assemble total likelihood function
                OBJmh = @(x)getObjective(x, ModelGuess, hasData, stateSpaces);
        
                %% Compute FIM for use in MH Proposal Function
                if strcmp(data,'simulated') %||testing
                    % Use the real FIM for the MH if we are using simulated data.
                    nCellsVec = zeros(1, controller.nInputs * controller.nT);
                    nCellsVec(1,((iInput-1)*controller.nT+1:iInput*controller.nT)) = ...
                        controller.nTotalCells(iInput,:);                    
                    FIM = ModelGuess{1}.evaluateExperiment(fimTrue, nCellsVec, ModelGuess{1}.fittingOptions.logPriorCovariance);
                else                                               
                    ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
                    fimResults = ModelGuess{iInput}.computeFIM([], fimScale);
                    
                    nCellsVec = controller.nTotalCells(iInput,:);
               
                    % Call function to assemble full FIM from cell
                    % counts and prior covariance information.
                    FIM = ModelGuess{1}.evaluateExperiment(...
                        fimResults, nCellsVec, ModelGuess{1}.fittingOptions.logPriorCovariance);
                end
            
                FIMfree = FIM{1}(ModelGuess{1}.fittingOptions.modelVarsToFit, ModelGuess{1}.fittingOptions.modelVarsToFit);
             
                if min(eig(FIMfree))<0.5
                    disp('Warning -- FIM has one or more small eigenvalues.  Reducing proposal in those directions. MH Convergence may be slow.')
                    FIMfree = FIMfree + 0.5*eye(length(FIMfree));
                end
             
                % Use FIM to define MHA proposal function (FIM -> covariance of MVN)
                covFree = FIMfree^-1;
                covFree = (covFree + covFree')/2;
                proposalDistribution = @(x) mvnrnd(x, covFree * FIMcovScale);
        
                if testing
                    MHFitOptions.progress = true;
                end
        
                %% Run MH Algorithm for tuning.
                [MHResults.mhSamples,MHResults.mhAcceptance,MHResults.mhValue,x0] = ...
                    ssit.parest.metropolisHastingsSample(log(newPars), ...
                    nSamplesMHtune, ...
                    'logpdf', OBJmh, 'proprnd', proposalDistribution, ...
                    'symmetric', true, ...
                    'thin', nThinMH, 'nchain', 1, 'burnin', nBurnMH, ...
                    'progress', MHFitOptions.progress, ...
                    'saveFileName', MHFitOptions.saveFile);
        
                delete(MHFitOptions.saveFile)
                
                mhTune = false;
                if max(abs(newPars-exp(x0))./newPars)>0  % If significant change in parameters   
                    % Convert best parameters to linear space.
                    newPars = exp(x0);
          
                    % Update models and statespaces
                    
                    ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
                    ModelGuess{iInput}.parameters(ModelGuess{iInput}.fittingOptions.modelVarsToFit, 2) = num2cell(newPars);
                    [fspSoln, ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
                    stateSpaces{iInput} = fspSoln.stateSpace;                       
                    
                    % Run another fminsearch starting at new parameter set.
                    % Call function to assemble total likelihood function
                    objFunc = @(x) -getObjective(x, ModelGuess, hasData, stateSpaces);
                    newPars = exp(fminsearch(objFunc, log(newPars), fitOptions));
        
                    % Update model parameters and statespace
                    
                    ModelGuess{iInput}.parameters(fitParameters, 2) = num2cell(newPars);
                    ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
                    [fspSoln, ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
                    stateSpaces{iInput} = fspSoln.stateSpace;
                                             
                    % Run another tuning round with new starting point.
                    mhTune = true;
        
                elseif MHResults.mhAcceptance<0.15
                    FIMcovScale = FIMcovScale/2;
        
                    % Run another tuning round with new starting point.
                    mhTune = true;
                end
            end % while mhTune ...
        
            disp('Starting New MH Chain')

            % Run new MH with tuned proposal.    
            [MHResults.mhSamples,MHResults.mhAcceptance,MHResults.mhValue,x0] = ...
                ssit.parest.metropolisHastingsSample(log(newPars), nSamplesMH,...
                'logpdf', OBJmh, 'proprnd', proposalDistribution, ...
                'symmetric', true, ...
                'thin', nThinMH, 'nchain', 1, 'burnin', nBurnMH, ...
                'progress', MHFitOptions.progress, ...
                'saveFileName', MHFitOptions.saveFile);
        
            % Update parameters if necessary.
            if max(abs(newPars-exp(x0))./newPars)>0  % If significant change in parameters
                % Convert best parameters to linear space.
                newPars = exp(x0);
        
                % Update models and statespaces
                
                ModelGuess{iInput}.fspOptions.fspTol = 1e-4;
                ModelGuess{iInput}.parameters(ModelGuess{iInput}.fittingOptions.modelVarsToFit,2) = num2cell(newPars);
                [fspSoln,ModelGuess{iInput}.fspOptions.bounds] = ModelGuess{iInput}.solve;
                stateSpaces{iInput} = fspSoln.stateSpace;                 
            end
        
            % Compute FIM for subsampling of MH results.
            J = floor(linspace(nSamplesMH/2,nSamplesMH,nFIMsamples));
            MHSamplesForFIM = exp(MHResults.mhSamples(J,:));
                      
            ModelGuess{iInput}.fspOptions.fspTol = 1e-8;
            fimResults = ModelGuess{iInput}.computeFIM([], fimScale, MHSamplesForFIM);                
     
            % FIM current experiment
            FIMCurrentExpt = ModelGuess{1}.totalFim(...
                fimResults, nCellsVec, ModelGuess{1}.fittingOptions.logPriorCovariance);
        
            switch lower(sampleType)
                case 'fimopt'
            % Find optimal NEXT experiment design given parameter sets
            nextExperiment = ModelGuess{1}.optimizeCellCounts(fimResults,numCellsPerExperiment,fimMetric,...
                [],nCellsVec,maxAvailable,'mean',...
                ModelGuess{1}.fittingOptions.logPriorCovariance,...
                incrementAdd);
            end
            FIMOptNextExpt = ModelGuess{1}.totalFim(fimResults,nextExperiment+nCellsVec,ModelGuess{1}.fittingOptions.logPriorCovariance);
        
            % Compute and Save Covariance from MH from CURRENT stage.
            parametersFound{iInput} = newPars;
            covLogMH{iInput} = cov(MHResults.mhSamples);
            covMH{iInput} = cov(exp(MHResults.mhSamples));
            MHResultsSaved{iInput} = MHResults;
            
            % Save FIM predictions for NEXT stage.
            FIMcurrentExptSaved{iInput} = FIMCurrentExpt;
            FIMcurrentExptTrueSaved{iInput} = FIMCurrentExpt_True;
            FIMpredNextExpt{iInput} = FIMOptNextExpt;
        
            % Save Predicted COV for NEXT stage.
            for jFIM=1:nFIMsamples
                covFIM_Prediction{iInput}{jFIM} = inv(FIMOptNextExpt{jFIM});
            end
        
            % Reshape next experiment design
            nextExperiment = reshape(nextExperiment',[nT,nInputs])';
            exptDesigns{iInput} = nextExperiment;
        
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
                f.Name = ['Current MH Results and Next FIM Prediction (Round ',num2str(iInput),')'];
                ModelGuess{1}.plotMHResults(MHResults,FIMOptNextExptReduced,fimScale,mhPlotScale,f)
        
                if iInput>1
                    f = figure;
                    f.Name = ['Current MH Results and Previous FIM Prediction (Round ',num2str(iInput),')'];
                    ModelGuess{1}.plotMHResults(MHResults,FIMpredNextExpt{iInput-1},fimScale,mhPlotScale,f)
                end
        
                f = figure;
                f.Name = ['Current MH Results and Perfect FIM Prediction (Round ',num2str(iInput),')'];
                ModelTrue{1}.plotMHResults(MHResults,FIMCurrentExpt_True,fimScale,mhPlotScale,f)
        
                figure;
                title(['Number of cells Measured at each Time Point (Round ',num2str(iInput),')']);
                bar(ModelGuess{1}.tSpan,nTotalCells,0.4, 'stacked')
                ylabel('Number of Cells Measured');
                xlabel('time [min]');
                ylim([0 300]);
            end
        
            % Save results
            save(saveExpName,'parametersFound','FIMcurrentExptSaved','FIMcurrentExptTrueSaved','covMH',...
                'covLogMH','exptDesigns','MHResultsSaved','FIMpredNextExpt')
        end
    end
end