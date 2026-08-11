classdef MetropolisHastingsAlgorithmRunner
    %MetropolisHastingsAlgorithmRunner tunes and runs the
    %Metropolis-Hastings algorithm for the purpose of exploring parameter
    %space and determining uncertainty in all parameter directions.

    properties (Access = private)
        IsTuning (1, 1) logical = true;
    end
    
    properties        
        GuessedModelsWithCombinedTimes (1, :) SSIT
        FIMCovarianceScale (1, 1) double {mustBePositive} = 1;
        FIMScale (1, 1) string = "log"
        FIMTrue
        FitOptions
        MHOptions (1, 1) MetropolisHastingsAlgorithmOptions       
        NumberOfSamplesForBurnIn (1, 1) uint64 {mustBePositive} = 1;
        NumberOfSamplesForProduction (1, 1) uint64 {mustBePositive} = 100;
        NumberOfSamplesForTuning (1, 1) uint64 {mustBePositive} = 100;
        NumberOfSamplesToThin (1, 1) uint64 {mustBePositive} = 2;
        ObjectiveFunction
        ObservationsMatrix (:, :) uint64 = []
        ParameterGuesses (1, :) double = [];
        ProposalDistribution
        StateSpaces (1, :) = [];
        TrueModelsWithCombinedTimes (1, :) SSIT
        UsingSimulatedData (1, 1) logical = true
    end

    methods (Access = private)
        function obj = computeFIMForUseInMH(obj)
            arguments
                obj (1, 1) MetropolisHastingsAlgorithmRunner
            end

            % Certain functions called here are effectively static, so we
            % will use the first guessed model as the caller.

            refModel = obj.GuessedModelsWithCombinedTimes(1);
            logPriorCovariance = ...
                refModel.fittingOptions.logPriorCovariance;

            numModels = length(obj.GuessedModelsWithCombinedTimes);
            timeSpan = MBExperimentDesigner.getAllTimes(...
                obj.GuessedModelsWithCombinedTimes);
            numTimes = length(timeSpan);

            nCellsVec = zeros(1, numModels * numTimes);

            if obj.UsingSimulatedData
                % Use the real FIM for the MH if we are using simulated 
                % data.
                
                for modelIdx = 1:numModels
                    range = MBExperimentDesigner.getArraySliceForModel(...
                        modelIdx, numTimes);
                    nCellsVec(1, range) = ...
                        obj.ObservationsMatrix(modelIdx, :);
                end
                FIM = refModel.evaluateExperiment(...
                    obj.FIMTrue, nCellsVec, logPriorCovariance);
            else
                fimResults = cell(numModels * numTimes, 1);
                for modelIdx = 1:numModels
                    curModel = obj.GuessedModelsWithCombinedTimes(...
                        modelIdx);
                    range = MBExperimentDesigner.getArraySliceForModel(...
                        modelIdx, numTimes);
                    fimResults(range, 1) = curModel.computeFIM(...
                        [], obj.FIMScale);
                    nCellsVec(1, range) = ...
                        obj.ObservationsMatrix(modelIdx, :);
                end

                % Call function to assemble full FIM from cell
                % counts and prior covariance information.
                FIM = refModel.evaluateExperiment(...
                    fimResults, nCellsVec, logPriorCovariance);
            end

            subsetVars = refModel.fittingOptions.modelVarsToFit;
            FIMfree = FIM{1}(subsetVars, subsetVars);

            if min(eig(FIMfree)) < 0.1
                disp("Warning -- FIM has one or more small " + ...
                    "eigenvalues.  Reducing proposal in those " + ...
                    "directions. MH Convergence may be slow.")
                FIMfree = FIMfree + 0.1 * eye(length(FIMfree));
            end

            % Use FIM to define MHA proposal function
            % (FIM -> covariance of MVN)
            covFree = FIMfree ^ -1;
            covFree = (covFree + covFree') / 2;
            obj.ProposalDistribution = ...
                @(x) mvnrnd(x, covFree * obj.FIMCovarianceScale);
        end % computeFIMForUseInMH

        function results = querySampler(obj)
            arguments
                obj (1, 1) MetropolisHastingsAlgorithmRunner
            end
        
            results = MetropolisHastingsAlgorithmResults();

            if obj.IsTuning
                numberOfSamples = obj.NumberOfSamplesForTuning;
            else
                numberOfSamples = obj.NumberOfSamplesForProduction;
            end
            
            [results.Samples, results.Acceptance, ...
                results.Value, results.Parameters] = ...
                ssit.parest.metropolisHastingsSample(...
                    log(obj.ParameterGuesses), ...
                    numberOfSamples, ...
                    'logpdf', ...
                    obj.ObjectiveFunction, ...
                    'proprnd', obj.ProposalDistribution, ...
                    'symmetric', true, ...
                    'thin', obj.NumberOfSamplesToThin, ...
                    'nchain', 1, ...
                    'burnin', obj.NumberOfSamplesForBurnIn, ...
                    'progress', obj.MHOptions.Progress, ...
                    'saveFileName', obj.MHOptions.SaveFile);

            if obj.IsTuning
                delete(obj.MHOptions.SaveFile)
            end
        end % querySampler

        function obj = reSolveAllModels(obj)
            arguments
                obj (1, 1) MetropolisHastingsAlgorithmRunner
            end

            oldModels = obj.GuessedModelsWithCombinedTimes;
            parameters = obj.ParameterGuesses;
            stateSpaces = [];

            parfor modelIdx = 1:length(oldModels)
                curModel = oldModels{modelIdx};

                curModel.fspOptions.fspTol = 1e-4;
                curModel.parameters(...
                    curModel.fittingOptions.modelVarsToFit, 2) = ...
                    num2cell(parameters);
                [curFspSoln, curModel.fspOptions.bounds] = ...
                    curModel.solve;
                stateSpaces{modelIdx} = curFspSoln.stateSpace;
            end

            obj.StateSpaces = stateSpaces;
        end % resolveAllModels

        function [obj, updateOccurred] = updateParameterGuesses(obj, ...
                proposedGuessesLogSpace, updateUnconditionally)
            arguments
                obj (1, 1) MetropolisHastingsAlgorithmRunner
                proposedGuessesLogSpace (1, :) double
                updateUnconditionally (1, 1) logical = false
            end

            newGuesses = exp(proposedGuessesLogSpace);

            if updateUnconditionally
                obj.ParameterGuesses = newGuesses;
                updateOccurred = true;
            else
                oldGuesses = obj.ParameterGuesses;                
                if max(abs(oldGuesses - newGuesses) ./ oldGuesses) > 0
                    obj.ParameterGuesses = newGuesses;
                    updateOccurred = true;
                end
            end

            if updateOccurred
                obj = obj.reSolveAllModels;
            end
        end % updateParameterGuesses    
    end % Private instance methods

    methods
        function [results, stateSpaces] = run(obj)
            arguments (Input)
                obj (1, 1) MetropolisHastingsAlgorithmRunner
            end
            arguments (Output)
                results (1, 1) MetropolisHastingsAlgorithmResults
                stateSpaces
            end

            %% Initialization of properties 

            if ~obj.UsingSimulatedData
                numModels = length(obj.GuessedModelsWithCombinedTimes);
                timeSpan = MBExperimentDesigner.getAllTimes(...
                    obj.GuessedModelsWithCombinedTimes);
                for modelIdx = 1:numModels
                    curModel = obj.GuessedModelsWithCombinedTimes(...
                        modelIdx);                    
                    curModel.tSpan = timeSpan;
                    curModel.fspOptions.fspTol = 1e-8;                    
                    obj.GuessedModelsWithCombinedTimes(modelIdx) = ...
                        curModel;
                end
            end
            %% MH Tuning then running.
        
            obj.IsTuning = true;
            while obj.IsTuning
        
                %% Metropolis Hastings
                disp('Starting New MH Chain for Tuning')                
                
                obj.MHOptions.Progress = false;
                obj.MHOptions.SaveFile = ...
                    ['TMPMHChain_',datFileName,'_',num2str(ind),'.mat'];
        
                % Delete old file if it exists.
                delete(MHFitOptions.saveFile)
        
                % Call function to assemble total likelihood function
                obj.ObjectiveFunction = ...
                    @(x) getLogLikelihoodOfDataGivenModels(...
                    x, obj.GuessedModelsWithCombinedTimes, ...
                    hasData, obj.StateSpaces);
        
                %% Compute FIM for use in MH Proposal Function
                obj = obj.computeFIMForUseInMH();
                       
                if testing
                    obj.MHOptions.Progress = true;
                end
        
                %% Run MH Algorithm for tuning.

                results = obj.querySampler();
               
                obj.IsTuning = false;

                % If a significant change occurred in the parameters, we
                % will continue tuning, after running another fminsearch
                % starting at the new parameter set.
                
                [obj, updateOccurred] = ...
                    obj.updateParameterGuesses(results.Parameters, false);
                if updateOccurred
                    obj.ObjectiveFunction = @(x) ...
                        -getLogLikelihoodOfDataGivenModels(...
                        x, obj.GuessedModelsWithCombinedTimes, ...
                        hasData, obj.StateSpaces);
                    newParameters = exp(fminsearch(...
                        obj.ObjectiveFunction, ...
                        obj.ParameterGuesses, obj.FitOptions));

                    obj = obj.updateParameterGuesses(newParameters, true);

                    % Run another tuning round with new starting point.
                    obj.IsTuning = true;
                elseif results.Acceptance < 0.15
                    obj.FIMCovarianceScale = obj.FIMCovarianceScale / 2;
        
                    % Run another tuning round with new starting point.
                    obj.IsTuning = true;
                end
            end % while obj.IsTuning
        
            %% Starting production run

            disp('Starting New MH Chain for Production')

            % Run new MH with tuned proposal, but first redefine the
            % objective function to take advantage of the most recently
            % calculated state spaces.

            obj.ObjectiveFunction = @(x) ...
                -getLogLikelihoodOfDataGivenModels(...
                x, obj.GuessedModelsWithCombinedTimes, ...
                hasData, obj.StateSpaces);
            results = obj.querySampler();

            % Update parameters if necessary.

            [obj, ~] = ...
                obj.updateParameterGuesses(results.Parameters, false);
            stateSpaces = obj.StateSpaces;
        end % run
    end % Public methods
end % class