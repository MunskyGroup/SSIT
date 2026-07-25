classdef MBExperimentDesigner
    %MBExperimentDesigner performs sequential model-based experiment design
    %   MBExperimentDesigner uses an SSIT model and empirical or simulated
    %   data to design and perform rounds of experimentation. In the case
    %   of empirical data, the user performs a round of experimentation by
    %   supplying a dataset that is novel to the designer; in the case of
    %   simulated data, the designer utilizes a ground-truth to simulate
    %   the outcome of a designed experiment.

    properties
        Configurations (1, :) MBExperimentConfiguration
        CumulativeNumbersOfObservations (1, :) uint64
        ExperimentRounds (1, :) MBExperimentRound
        FIMScale (1, 1) string = "log"
        GroundTruthModel (1, 1) SSIT
        GuessedModel (1, 1) SSIT
        MaxFitIterations (1, 1) uint64 {mustBePositive} = 1;
        MHOptions (1, 1) MetropolisHastingsAlgorithmOptions
        NextExperimentDesign (1, 1) MBExperimentDesign
        NumberOfFIMSamples (1, 1) uint64 {mustBePositive} = 1
        NumberOfFitRounds (1, 1) uint64 {mustBePositive} = 1
        NumberOfMHSamplesForBurnIn (1, 1) uint64 {mustBePositive} = 1;
        NumberOfMHSamplesForProduction (1, 1) uint64 {mustBePositive} = 100;
        NumberOfMHSamplesForTuning (1, 1) uint64 {mustBePositive} = 100;
        NumberOfMHSamplesToThin (1, 1) uint64 {mustBePositive} = 2;
        Property1
        Strategy (1, 1) MBAbstractExperimentDesignStrategy = ...
            MBExperimentRandomDesignStrategy();        
        UseEmpiricalData (1, 1) logical = false
    end

    properties (Dependent)
        FitOptions
        GroundTruthModelsWithCombinedTimes (1, :) SSIT
        GroundTruthModelsWithIndividualTimes (1, :) SSIT
        GuessedModelsWithCombinedTimes (1, :) SSIT
        GuessedModelsWithIndividualTimes (1, :) SSIT        
    end

    methods (Static)
        function timeSpan = getAllTimes(configurations)
            arguments
                configurations (1, :) MBExperimentConfiguration
            end

            timeSpan = [];
            for configIdx = 1:length(configurations)
                timeSpan = [timeSpan ...
                    configurations(configIdx).TimeConfigurable.Values];                
            end
            timeSpan = unique(sort(timeSpan));
        end

        function range = getArraySliceForModel(modelIdx, numTimes)
            arguments
                modelIdx (1, 1) int64 {mustBePositive}
                numTimes (1, 1) int64 {mustBePositive}
            end

            range = ((modelIdx - 1) * numTimes + 1) : (modelIdx * numTimes);
        end
    end % Public static methods

    % TO DO: collectSummaryOfDesignedRound needs 
    % fimResults, nextExperiment, nCellsVec,
    % FIMOptNextExpt

    methods (Access = private)
        function round = collectSummaryOfDesignedRound(obj, round)
            arguments
                obj (1, 1) MBExperimentDesigner
                round (1, 1) MBExperimentRound
            end
            
            round.FIMPredNextExpt = totalFim(...
                round.FIMResults, ...
                nextExperiment+obj.CumulativeNumbersOfObservations, ...
                obj.GuessedModel.fittingOptions.logPriorCovariance);

            % Compute and Save Covariance from MH from CURRENT stage.
            
            round.ParametersFound = round.MHResults.Parameters;
            round.CovarianceLogMH = cov(round.MHResults.Samples);
            round.CovarianceMH = cov(exp(round.MHResults.Samples));            

            % Save Predicted COV for NEXT stage.
            covarianceFIM_Prediction = ...
                createArray(1, obj.NumberOfFIMSamples, "cell");
                       
            for idxFIMSample = 1 : obj.NumberOfFIMSamples
                covarianceFIM_Prediction{idxFIMSample} = ...
                    inv(round.FIMOptNextExpt{idxFIMSample});
            end
            round.CovarianceFIM_Prediction = covarianceFIM_Prediction;
        end % collectSummaryOfDesignedRound

        function round = computeFIMsForDesignedRound(obj, round)
            arguments
                obj (1, 1) MBExperimentDesigner
                round (1, 1) MBExperimentRound
            end
           
            % Compute FIM for subsampling of MH results; choose the later
            % half of the production run.
          
            J = floor(linspace(...
                obj.NumberOfMHSamplesForProduction / 2, ...
                obj.NumberOfMHSamplesForProduction, ...
                obj.NumberOfFIMSamples));
            MHSamplesForFIM = exp(round.MHResults.Samples(J, :));

            models = obj.GuessedModelsWithCombinedTimes;
            numModels = length(models);
            timeSpan = getAllTimes(obj.GuessedModelsWithCombinedTimes);
            numTimes = length(timeSpan);

            % It is simpler and not excessively more costly to compute the
            % FIM for all times for all models, even if not all
            % combinations of models and times are experimentally
            % available. We therefore create a copy of each model that
            % includes all times across the design space, for the purpose
            % of FSP solution so that the appropriate FIMs can be
            % calculated. However, we do not need to store these modified
            % models.

            fimResults = cell(...
                numModels * numTimes, obj.NumberOfFIMSamples);
            for modelIdx = 1:numModels
                curModel = obj.GuessedModelsWithCombinedTimes(modelIdx);

                curModel.tSpan = timeSpan; % Include all times
                curModel.fspOptions.fspTol = 1e-8;

                range = MBExperimentDesigner.getArraySliceForModel(...
                    modelIdx, numTimes);

                fimResults(range, :) = curModel.computeFIM(...
                    [], obj.FIMScale, MHSamplesForFIM);                
            end
            round.FIMResults = fimResults;
          
            % Compute and save FIM predictions for NEXT stage.

            % FIM current experiment
            round.FIMCurrentExpt = totalFim(fimResults, ...
                obj.CumulativeNumbersOfObservations, ...
                obj.GuessedModel.fittingOptions.logPriorCovariance);

            % True FIM for current experiment
            round.FIMCurrentExpt_True = totalFim(fimTrue, ...
                obj.CumulativeNumbersOfObservations, ...
                obj.GuessedModel.fittingOptions.logPriorCovariance);
        end

        function round = designNextExperiment(obj, round)
            arguments
                obj (1, 1) MBExperimentDesigner
                round (1, 1) MBExperimentRound
            end

            [round.NextExperimentDesign, ...
                cellVecForOptimalFIMCalculation] = ...
                ...
                obj.Strategy.apportionObservations(...
                    round.NextExperimentDesign);

            round.FIMOptNextExpt = totalFim(round.FIMResults, ...
                cellVecForOptimalFIMCalculation + ...
                obj.CumulativeNumbersOfObservations, ...
                obj.GuessedModel.fittingOptions.logPriorCovariance);            
        end

        function [obj, round] = setupAndRunMetropolisHastings(obj, round)
            arguments
                obj (1, 1) MBExperimentDesigner
                round (1, 1) MBExperimentRound
            end
            
            runner = MetropolisHastingsAlgorithmRunner();
            
            runner.FIMScale = obj.FIMScale;

            % FIMTrue               

            runner.FitOptions = obj.FitOptions;       
            
            runner.GuessedModelsWithCombinedTimes = ...
                obj.GuessedModelsWithCombinedTimes;

            runner.MHOptions = obj.MHOptions;

            runner.NumberOfSamplesForBurnIn = ...
                obj.NumberOfMHSamplesForBurnIn;
            runner.NumberOfSamplesForProduction = ...
                obj.NumberOfMHSamplesForProduction;
            runner.NumberOfSamplesForTuning = ...
                obj.NumberOfMHSamplesForTuning;
            runner.NumberOfSamplesToThin = ...
                obj.NumberOfMHSamplesToThin;

            runner.ParameterGuesses = [obj.GuessedModel.parameters{...
                obj.GuessedModel.fittingOptions.modelVarsToFit, 2}];
            
            % StateSpaces (1, :) = [];         

            runner.TrueModelsWithCombinedTimes = ...
                obj.GroundTruthModelsWithCombinedTimes;
            
            runner.UsingSimulatedData = ~obj.UseEmpiricalData;
            
            % Although the runner also provides up-to-date state spaces for
            % models, we do not utilize them downstream from here. We will,
            % however, of course update the parameters in the guessed model
            % based on the M-H results.

            [round.MHResults, ~] = runner.run();
            
            obj.GuessedModel.parameters(...
                obj.GuessedModel.fittingOptions.modelVarsToFit, 2) = ...
                num2cell(round.MHResults.Parameters);            
        end
    end % Private methods

    methods
        function obj = MBExperimentDesigner(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function [nextDesign, obj] = designNextRound(obj)
            %designNextRound determines an experiment design, i.e., an
            %apportionment of available observations across candidate
            %experiment configurations, accordingly to the design strategy
            %associated with the designer.

            round = MBExperimentRound();
            round.CumulativeNumbersOfObservations = ...
                obj.CumulativeNumbersOfObservations;
            
            % 1. Fit new parameter values to updated data
            % 2. Tune and run Metropolis-Hastings; update parameters to
            % those of MLE estimator

            % 3. Compute FIMs from Metropolis-Hastings results
            
            round = obj.computeFIMsForDesignedRound(round);
            
            % 4. Design the next experiment, i.e., apportion available 
            % observations to configurations.

            obj.Strategy = obj.Strategy.incorporateRoundDetails(round);
            round = obj.designNextExperiment(round);
            
            % 5. Collect remaining statistics and summary of next round

            round = obj.collectSummaryOfDesignedRound(round);

            % 6. Finalize by including the round in the sequence of rounds
            % and saving the experiment design for ready use in the next
            % invocation of performNextRound.            

            obj.ExperimentRounds(end + 1) = round;

            obj.NextExperimentDesign = round.NextExperimentDesign;
            nextDesign = obj.NextExperimentDesign;
        end % designNextRound

        function options = get.FitOptions(obj)
            options = optimset('Display', 'iter', ...
                'MaxIter', obj.MaxFitIterations);
        end

        function [nextRoundResults, obj] = performNextRound(obj)
        end % performNextRound

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end % Public methods
end