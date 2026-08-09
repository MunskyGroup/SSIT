classdef MBExperimentDesigner < handle
    %MBExperimentDesigner performs sequential model-based experiment design
    %   MBExperimentDesigner uses a model and empirical or simulated data
    %   to design and perform rounds of experimentation. In the case of
    %   empirical data, the user performs a round of experimentation by
    %   supplying a dataset that is novel to the designer; in the case of
    %   simulated data, the designer utilizes a ground-truth model to 
    %   simulate the outcome of a designed experiment.

    properties (Constant)
        TABLE_ROUND_COLUMN = "UTILIZED_IN_PERFORMED_ROUND"
    end

    properties (SetAccess = private)
        CacheOfGroundTruthModelsWithCombinedTimes (1, 1) ...
            MBExperimentDependentModelsCache
        CacheOfGroundTruthModelsWithIndividualTimes (1, 1) ...
            MBExperimentDependentModelsCache
        CacheOfGuessedModelsWithCombinedTimes (1, 1) ...
            MBExperimentDependentModelsCache
        CacheOfGuessedModelsWithIndividualTimes (1, 1) ...
            MBExperimentDependentModelsCache
        CumulativeNumbersOfObservations (1, :) uint64
        DesignedExperimentRounds (1, :) MBExperimentRound = ...
            createArray(0, 0, "MBExperimentRound")
        EmpiricalDataFiles (1, :) string = createArray(0, 0, "string")       
        EmpiricalDataTable table = table()
        GlobalRNGStream RandStream = []
        GroundTruthModel (1, 1) MBExperimentModel = []
        LocalRNGStream RandStream = []
        PerformingRoundNumber (1, 1) uint64 = 0        
        ReadyToPerformNextRound (1, 1) logical = false
        SimulatedDataFiles (1, :) string = createArray(0, 0, "string")
        SimulatedDataTable table = table()
    end % Read-only properties

    properties
        % Note that the class of Configurations property is not restricted 
        % to MBExperimentConfiguration here; the setter accepts multiple
        % classes of arguments and handles necessary conversions.
        Configurations (1, :) = ...
            createArray(0, 0, "MBExperimentConfiguration")
        ErrorOnInsufficientAvailableObservations (1, 1) logical = false       
        FIMScale (1, 1) string = "log"       
        GuessedModel (1, 1) MBExperimentModel = []
        MaxFitIterations (1, 1) uint64 {mustBePositive} = 1
        MHOptions (1, 1) MetropolisHastingsAlgorithmOptions
        ModelToDataColumnsMap (1, 1) dictionary = ...
            configureDictionary("string", "string")
        NextExperimentDesign (1, 1) MBExperimentDesign
        NumberOfFIMSamples (1, 1) uint64 {mustBePositive} = 1
        NumberOfFitRounds (1, 1) uint64 {mustBePositive} = 1
        NumberOfMHSamplesForBurnIn (1, 1) uint64 {mustBePositive} = 1
        NumberOfMHSamplesForProduction (1, 1) uint64 {mustBePositive} = 100
        NumberOfMHSamplesForTuning (1, 1) uint64 {mustBePositive} = 100
        NumberOfMHSamplesToThin (1, 1) uint64 {mustBePositive} = 2
        Strategy (1, 1) MBAbstractExperimentDesignStrategy = ...
            MBExperimentRandomDesignStrategy()        
        UseEmpiricalData (1, 1) logical = true
    end

    properties (Dependent)
        FitOptions
        GroundTruthModelsWithCombinedTimes (1, :) MBExperimentModel
        GroundTruthModelsWithIndividualTimes (1, :) MBExperimentModel
        GuessedModelsWithCombinedTimes (1, :) MBExperimentModel
        GuessedModelsWithIndividualTimes (1, :) MBExperimentModel                
    end

    methods (Static, Access = private)
        function configurations = resolveConfigurations(configs)
            % This helper method "resolves" the input argument into a list
            % of configurations:
            % If a list of configurations proper (i.e., 
            % MBExperimentConfiguration) is provided, the list is cloned.
            % If a list of CONFIGURABLES (i.e., 
            % MBAbstractExperimentConfigurable) is provided, the 
            % configurables will be "multiplied" into distinct 
            % configurations, each with a single value for each
            % configurable.

            arguments
                configs (1, :)
            end

            if isa(configs, "MBExperimentConfiguration")
                configurations = configs;
            elseif isa(configs, "MBAbstractExperimentConfigurable")
                configurations = ...
                    MBExperimentDesigner.multiplyConfigurables(configs);
            else
                error("Attempt to create or modify experiment " + ...
                    "configurations using an unsupported class: " + ...
                    class(configs));
            end
        end % resolveConfigurations
    end % Private static methods

    methods (Static)
        function configsOut = combineConfigurationsByTime(configsIn)
            %combineConfigurationsByTime accepts a list of experiment
            %configurations and outputs an equivalent list of experiment
            %configurations such that
            %   1. Each output configuration contains a unique non-time
            %   configuration (set of values for input expressions and 
            %   designed parameters).
            %   2. As a result of the above, all time configurables
            %   corresponding to each non-time configuration are combined
            %   into a single time configurable.            
            arguments
                configsIn (1, :) MBExperimentConfiguration
            end

            d = configureDictionary(...
                "MBExperimentNonTimeConfiguration", ...
                "MBExperimentTimeConfigurable");
            
            for inIdx = 1:length(configsIn)
                ntc = configsIn(inIdx).NonTimeConfiguration;
                tc = configsIn(inIdx).TimeConfigurable;
                if d.isKey(ntc)
                    % There is already a corresponding time configurable
                    % for this non-time configuration. Append that existing
                    % configurable to the incoming one and update the value
                    % accordingly.

                    d(ntc) = d(ntc) + tc;
                else
                    % There is no existing corresponding time configurable
                    % for this non-time configuration, so introduce a new
                    % key-value pair.

                    d(ntc) = tc;
                end
            end % combineConfigurationsByTime

            % Having now generated a dictionary mapping all non-time
            % configurations to combined time configurables, generate the
            % output configurations by pairing the keys and values:

            configsOut = createArray(1, d.numEntries(), ...
                "MBExperimentConfiguration");
            keys = d.keys();
            for outIdx = 1:length(keys)
                newConfig = MBExperimentConfiguration();
                ntc = keys(outIdx);
                newConfig.NonTimeConfiguration = ntc;
                newConfig.TimeConfigurable = d(ntc);
                configsOut(outIdx) = newConfig;               
            end                          
        end

        function timeSpan = getAllTimes(configurations)
            arguments
                configurations (1, :) MBExperimentConfiguration
            end

            if ~isempty(configurations)
                % Sum all the configurations' TimeConfigurable; this has
                % the effect of combining all their times into the Value of
                % the sum.

                tc = configurations(1).TimeConfigurable;
                for configIdx = 2:length(configurations)
                    tc = tc + configurations(configIdx).TimeConfigurable;
                end
                timeSpan = tc.Values;
            else
                timeSpan = []; % No configurations means no times...
            end
        end % getAllTimes

        function range = getArraySliceForModel(modelIdx, numTimes)
            arguments
                modelIdx (1, 1) int64 {mustBePositive}
                numTimes (1, 1) int64 {mustBePositive}
            end

            range = ((modelIdx - 1) * numTimes + 1) : (modelIdx * numTimes);
        end

        function configurations = multiplyConfigurables(configurables)
            arguments
                configurables (1, :) MBAbstractExperimentConfigurable
            end

            if isempty(configurables)
                configurations = createArray(0, 0, ...
                    "MBExperimentConfiguration");
            else                
                if isscalar(configurables)
                    % We multiply the single configurable into N
                    % configurations, where N is the number of possible
                    % values of the configurable:

                    N = configurables.numberOfValues();

                    configurations = createArray(...
                        1, N, "MBExperimentConfiguration");
                    
                    multipliedConfigurables = configurables.multiply();

                    for configIdx = 1:N
                        configurations(configIdx) = ...
                            configurations(configIdx).setConfigurable(...
                                multipliedConfigurables(configIdx));
                    end                    
                else
                    % We operate recursively by first multiplying all 
                    % configurables except the last one; suppose that this 
                    % results in M configurations, whereas the last 
                    % configurable has N possible values. Then we will have
                    % M*N configurations after multiplying by the last 
                    % configurable.

                    preMultipliedConfigurations = ...
                        MBExperimentDesigner.multiplyConfigurables(...
                            configurables(1:end-1));          
                    
                    M = length(preMultipliedConfigurations);

                    lastConfigurables = configurables(end).multiply();
                    N = length(lastConfigurables);           

                    configurations = ...
                        createArray(1, M * N, "MBExperimentConfiguration");

                    % We will "stripe" the configurations such that the
                    % first N correspond to the first "pre-multiplied"
                    % configuration with the incorporation of the N values
                    % of the "last" configurable, the second N correspond
                    % to the second "pre-multiplied" configuration with the
                    % incorporation of the "last" configurable, etc.

                    configurationIdx = 1;
                    for preMultipliedConfigIdx = 1:M
                        preConfig = preMultipliedConfigurations(...
                            preMultipliedConfigIdx);

                        for lastConfigurableIdx = 1:N
                            configurable = lastConfigurables(...
                                lastConfigurableIdx);

                            configurations(configurationIdx) = ...
                                preConfig.setConfigurable(configurable);

                            configurationIdx = configurationIdx + 1;
                        end
                    end % [penultimate configurables]
                end % [configurables contains 2+]
            end % [configurables is non-empty]
        end % multiplyConfigurables        
    end % Public static methods

    methods (Access = private)
        function buildDataTable(obj, dataFiles)
            arguments
                obj (1, 1) MBExperimentDesigner
                dataFiles (1, :) string
            end

            if obj.UseEmpiricalData
                existingDataFiles = obj.EmpiricalDataFiles;
            else
                existingDataFiles = obj.SimulatedDataFiles;
            end

            % Subset the list of incoming data files to those not already
            % seen.

            dataFileIndicesNotAlreadySeen = ...
                createArray(1, length(dataFiles), "logical");
            for fileIdx = 1:length(dataFiles)
                dataFileIndicesNotAlreadySeen(fileIdx) = ~any(...
                    strcmp(existingDataFiles, dataFiles(fileIdx)));
            end
            dataFiles = dataFiles(dataFileIndicesNotAlreadySeen);

            % Create a table from all new incoming data files, if there are
            % any. Then create a new column indicating the round, if any,
            % in which the observation (row) was utilized. This is
            % initialized to NaN rather than to zero for clarity and for
            % fast lookup of which observations have (or haven't) been
            % utilized.

            if ~isempty(dataFiles)
                ds = tabularTextDatastore(dataFiles);
                T = ds.readall;
                T.(obj.TABLE_ROUND_COLUMN) = NaN(height(T), 1);

                if obj.UseEmpiricalData
                    % Use the ModelToDataColumnsMap to rename columns to
                    % match corresponding features of the model
                    % (parameters, species, input expressions).

                    varNames = string(T.Properties.VariableNames);
                    d = obj.ModelToDataColumnsMap;
                    k = d.keys;
                    for keyIdx = 1:d.numEntries()
                        modelColumn = k(keyIdx);
                        dataColumn = d(modelColumn);                        
                        varNames(strcmp(varNames, dataColumn)) = ...
                            modelColumn;
                    end
                    T.Properties.VariableNames = varNames;
                end

                % Append the new table to the corresponding existing table
                % and the new data files to the corresponding list of
                % existing data files.

                if obj.UseEmpiricalData
                    obj.EmpiricalDataFiles = [existingDataFiles dataFiles];
                    obj.EmpiricalDataTable = [obj.EmpiricalDataTable; T];
                else
                    obj.SimulatedDataFiles = [existingDataFiles dataFiles];
                    obj.SimulatedDataTable = [obj.SimulatedDataTable; T];
                end
            end % [New incoming data files exist]
        end % buildDataTable


        % TO DO: collectSummaryOfDesignedRound needs 
        % nextExperiment, nCellsVec   

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
        end % computeFIMsForDesignedRound

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
        end % designNextExperiment

        function [models, configs] = multiplyModel(obj, ...
                model, configs, combineTimes, refitToNewData, cache)
            %multiplyModel accepts a source model, a list of experiment
            %configurations, and a Boolean indicating whether measurement
            %times should be combined for otherwise similar configurations,
            %and it generates a cache of models, configurations, and
            %related properties such that
            %   1. Each output configuration is the result of combining the
            %   measurement times in the corresponding source
            %   configurations (if indicated),
            %   2. Each output model is the result of
            %       a. Applying the corresponding configuration to the
            %       source model, including the measurement times;
            %       b. Loading all existing relevant data into the model:
            %           i. Data must match the configuration, and
            %           ii. Data must have been utilized in a performed
            %           experiment round.
            %       c. Fitting the collection of output models to the data
            %       and updating all models' fitted parameters accordingly,
            %       IF the caller has indicated that models should be refit
            %       (which should NOT be done for ground-truth models).

            arguments                
                obj (1, 1) MBExperimentDesigner
                model (1, 1) MBExperimentModel
                configs (1, :) MBExperimentConfiguration
                combineTimes (1, 1) logical
                refitToNewData (1, 1) logical
                cache (1, 1) MBExperimentDependentModelsCache = [];
            end

            if ~isempty(cache)
                % If an input cache has been provided, assume nothing needs
                % to be regenerated/reloaded/refit, then begin to
                % check:
                regenerateConfigs = false;
                regenerateModels = false;
                reloadData = false;

                % If the configs have changed, then everything will need to
                % be regenerated:
                %   1. Models are generated by applying each config to the
                %   source model.
                %   2. Data are loaded according to the config.
                %   3. Models are fit according to the loaded data.

                if length(configs) ~= length(cache.Configs)
                    regenerateConfigs = true;
                else
                    for configIdx = 1:length(configs)
                        if ~(configs(configIdx) == ...
                                cache.Configs(configIdx))
                            regenerateConfigs = true;
                            break
                        end
                    end
                end % [Regenerate configs?]

                if regenerateConfigs
                    regenerateModels = true;
                    reloadData = true;
                end

                % Even if the configs haven't changed, models need to be
                % regenerated if the source model has changed:

                if ~(model == cache.SourceModel)
                    regenerateModels = true;
                end

                % Even if neither the models nor configs have changed, data
                % need to be reloaded if 
                %   1. A different type of data is being used in this round
                %   (simulated versus empirical), and therefore a different
                %   source data table,
                %   OR
                %   2. More data have been acquired (as indicated by an
                %   outdated number of performed experiment round).

                if obj.UseEmpiricalData ~= cache.UseEmpiricalData
                    reloadData = true;
                elseif obj.PerformingRoundNumber > ...
                        cache.LastPerformedExperimentRound
                    reloadData = true;
                end               
            else % [Input cache provided]
                % If no input cache is provided, everything needs to be
                % regenerated/reloaded/recalculated.

                cache = MBExperimentDependentModelsCache();
                regenerateConfigs = true;
                regenerateModels = true;
                reloadData = true;
            end % [No input cache provided]

            % If data need to be reloaded, then models that should be
            % refit to new data will need to be refit accordingly: 

            refitModels = refitToNewData && reloadData; 

            % Regenerate configs if necessary:

            if regenerateConfigs
                if combineTimes
                    configs = combineConfigurationsByTime(configs);
                end
            else
                configs = cache.Configs;
            end
            numConfigs = length(configs); % Also the number of models
            
            % Regenerate models if necessary:

            if regenerateModels
                models = createArray(1, numConfigs, "MBExperimentModel");

                for configIdx = 1:numConfigs
                    models(configIdx) = ...
                        configs(configIdx).applyToModel(model);
                end

                haveData = zeros(1, numConfigs, "logical");
                stateSpaces = cell(1, numConfigs);
            else
                models = cache.Models;
                haveData = cache.ModelsHaveData;
                stateSpaces = cache.ModelStateSpaces;
            end

            % Reload data if necessary:

            if reloadData
                % Get the correct data to (re)load:

                if obj.UseEmpiricalData
                    data = obj.EmpiricalDataTable;
                else
                    data = obj.SimulatedDataTable;
                end

                mapToUse = obj.ModelToDataColumnsMap;

                % In parallel, subset the data for each model, solve the
                % model using FSP, and assign the solved model (including
                % its solution, which contains the state space) back to the
                % list of models.
               
                parfor modelIdx = 1:length(models)
                    curConfig = configs(modelIdx);
                    curData = curConfig.applyToData(data);
                    haveData(modelIdx) = height(curData) > 0;
                    if haveData(modelIdx)
                        curModel = models(modelIdx);

                        curModel = curModel.loadData(curData, mapToUse);
                        curModel.fspOptions.fspTol = 1e-4;
                        curModel = curModel.solve(solver = "fsp");

                        models(modelIdx) = curModel;
                        stateSpaces(modelIdx) = ...
                            curModel.Solutions.stateSpace;
                    end                 
                end
            end % [Reload data]

            if refitModels
                fitParameters = model.fittingOptions.modelVarsToFit;
                parameters = model.parameters(fitParameters, 2);

                for fitRoundIdx = 1:obj.NumberOfFitRounds

                    % Find the best parameter values, i.e., those that
                    % maximize the (log-)likelihood of all data given all
                    % models:

                    objective = @(x) -getLogLikelihoodOfDataGivenModels(...
                        x, models, haveData, stateSpaces);                    
                    parameters = exp(fminsearch(...
                        objective, log(parameters), obj.FitOptions));

                    % Update and resolve all models:

                    model.parameters(fitParameters, 2) = ...
                        num2cell(parameters);

                    parfor modelIdx = 1:numConfigs
                        if hasData(modelIdx)
                            curModel = model(modelIdx);

                            curModel.parameters(fitParameters, 2) = ...
                                num2cell(parameters);
                            curModel.fspOptions.fspTol = 1e-8;

                            curModel = curModel.solve(solver = "fsp");

                            models(modelIdx) = curModel;
                            stateSpaces(modelIdx) = ...
                                curModel.Solutions.stateSpace;
                        end                       
                    end % [Model]
                end % [Fit round]
            end % [Refit models]

            % Having completed the multiplication, update the cache:

            cache.Configs = configs;
            cache.LastPerformedExperimentRound = obj.PerformingRoundNumber;
            cache.Models = models;
            cache.ModelsHaveData = haveData;
            cache.ModelStateSpaces = stateSpaces;
            cache.SourceModel = model;
            cache.UseEmpiricalData = obj.UseEmpiricalData;        
        end % multiplyModel
      
        function round = setupAndRunMetropolisHastings(obj, round)
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
        end % setupAndRunMetropolisHastings

        function sampleObservationsPerDesign(obj)
            arguments
                obj (1, 1) MBExperimentDesigner
            end

            if obj.UseEmpiricalData
                T = obj.EmpiricalDataTable;
            else
                T = obj.SimulatedDataTable;
            end

            unsampledRows = isnan(T.(obj.TABLE_ROUND_COLUMN));

            design = obj.NextExperimentDesign;
            configs = design.Configurations;

            numConfigs = length(configs);
            for configIdx = 1:numConfigs
                curConfig = configs(configIdx);
                curNumObservations = design.getObservationsForConfiguration(...
                    curConfig);
                if curNumObservations > 0
                    % If the design involves making observations for this
                    % configuration, then find
                    %   1. The subset of the data pertaining to this
                    %   configuration.
                    %   2. The UNSAMPLED subset of the first subset (by
                    %   logical conjunction).
                    %   3. A third subset of the second subset equal in
                    %   length to the number of desired observations.

                    rowsForCurConfig = curConfig.findSubsetOfData(T);
                    unsampledRowsForCurConfig = ...
                        unsampledRows & rowsForCurConfig;
                    rowsToSample = find(...
                        unsampledRowsForCurConfig, curNumObservations);

                    if obj.ErrorOnInsufficientAvailableObservations && ...
                            length(rowsToSample) < curNumObservations
                        error("Desired " + curNumObservations + ...
                            " but only " + length(rowsToSample) + ...
                            " were available for configuration " + ...
                            curConfig.FilenameString)
                    end

                    % Having found the desired number of appropriate rows
                    % to sample, sample them by assigning their values in
                    % the "round" column to the number of the round being
                    % performed.

                    T.(obj.TABLE_ROUND_COLUMN)(rowsToSample) = ...
                        obj.PerformingRoundNumber;
                end % [observations desired for this configuration]
            end % [for each configuration]

            % Finally, having sampled sufficiently many desired and
            % available observations for each configuration, assign the
            % updated table back to the appropriate property:

            if obj.UseEmpiricalData
                obj.EmpiricalDataTable = T;
            else
                obj.SimulatedDataTable = T;
            end
        end % sampleObservationsPerDesign

        function outputFiles = simulateSystem(obj)
            arguments
                obj (1, 1) MBExperimentDesigner
            end

            % From the design, obtain a map from the non-time
            % configurations to the maximum number of observations at any
            % corresponding time point. We do this to ensure that, because
            % each simulated model combines all time points for a given
            % non-time configuration, we will have sufficiently many
            % observations for each configuration (non-time plus time).

            numObservationsMap = ...
                getMostObservationsAtAnyTimeForNonTimeConfigurations(...
                obj.NextExperimentDesign);

            % Obtain a list of true models, one for experiment
            % configuration after having combined measurement times. Note
            % that there is a 1-1 map between models and configs.

            % The first method handles updating (if necessary) the cache 
            % and returning the models. We can then query the cache for the
            % configs:

            models = obj.GroundTruthModelsWithCombinedTimes;
            configs = ...
                obj.CacheOfGroundTruthModelsWithCombinedTimes.Configs;

            numConfigs = length(configs);
            maxObservations = zeros(1, numConfigs, "uint64");
            for configIdx = 1:numConfigs
                curNonTimeConfig = configs(configIdx).NonTimeConfiguration;

                % If the non-time configuration is not a key of the map,
                % that means that design specifies no observations at any
                % time for that configuration, so the value should remain
                % zero.

                if numObservationsMap.isKey(curNonTimeConfig)
                    maxObservations(configIdx) = ...
                        numObservationsMap(curNonTimeConfig);
                end
            end

            numModels = length(models);
            outputFiles = createArray(1, numModels, "string");
            for modelIdx = 1:numModels
                outputFiles(modelIdx) = "SimulatedData_" + ...
                    "Round_" + obj.PerformingRoundNumber + ...
                    "Model_" + modelIdx + ".csv";
            end

            parfor modelIdx = 1:numModels
                curModel = models(modelIdx);                
                curModel.ssaOptions.nSimsPerExpt = ...
                    maxObservations(modelIdx);
                curModel.ssaOptions.Nexp = 1;

                if curModel.ssaOptions.nSimsPerExpt ~= 0
                    curModel.sampleDataFromFSP(...
                        [], ... % Leave empty to use existing FSP solution
                        outputFiles(modelIdx));
                end
            end % [for each model]
        end % simulateSystem
    end % Private methods

    methods
        function excludeConfigurations(obj, configs)
            arguments
                obj (1, 1) MBExperimentDesigner
                configs (1, :)
            end

            configsToExclude = resolveConfigurations(configs);
            for excludeConfigIdx = 1:length(configsToExclude)
                configToExclude = configsToExclude(excludeConfigIdx);

                % Remove the config from the list of configurations:

                for configIdx = length(obj.Configurations):-1:1
                    if obj.Configurations(configIdx) == configToExclude
                        % The Configurations setter validates that the
                        % configurations are all unique, so at most one
                        % configuration should match each configuration to
                        % exclude:

                        obj.Configurations(configIdx) = [];
                        break
                    end
                end

                % Remove the config from the next experiment design:

                obj.NextExperimentDesign.excludeConfiguration(...
                    configToExclude);
            end
        end % excludeConfigurations

        function nextDesign = designNextRound(obj)
            %designNextRound determines an experiment design, i.e., an
            %apportionment of available observations across candidate
            %experiment configurations, accordingly to the design strategy
            %associated with the designer.

            round = MBExperimentRound();
            round.CumulativeNumbersOfObservations = ...
                obj.CumulativeNumbersOfObservations;

            % First, switch to using our local RNG stream, if available.
            % This is the last point before randomness will be utilized,
            % e.g., in sampling parameter space. It is also the last point
            % where we are guaranteed not to error.

            if ~isempty(obj.LocalRNGStream)
                obj.GlobalRNGStream = ...
                    RandStream.setGlobalStream(obj.LocalRNGStream);
            else
                obj.GlobalRNGStream = RandStream.getGlobalStream();
            end

            try
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

                % 6. Finalize by including the round in the sequence of 
                % rounds and saving the experiment design for ready use in
                % the next invocation of performNextRound.            

                obj.DesignedExperimentRounds(end + 1) = round;

                obj.NextExperimentDesign = round.NextExperimentDesign;
                nextDesign = obj.NextExperimentDesign;

                obj.ReadyToPerformNextRound = true;

                RandStream.setGlobalStream(obj.GlobalRNGStream);
            catch ME
                % Reset the global stream and rethrow the error:
                RandStream.setGlobalStream(obj.GlobalRNGStream);
                rethrow ME
            end % [Design, with potential faults]                  
        end % designNextRound

        function options = get.FitOptions(obj)
            options = optimset('Display', 'iter', ...
                'MaxIter', obj.MaxFitIterations);
        end

        function models = get.GroundTruthModelsWithCombinedTimes(obj)
            obj.CacheOfGroundTruthModelsWithCombinedTimes = ...
                obj.multiplyModel(obj.GroundTruthModel, ...
                    obj.Configurations, ...
                    true, ... % Combine times
                    false, ... % Ground-truth models, so DON'T refit 
                    obj.CacheOfGroundTruthModelsWithCombinedTimes);
            models = ...
                obj.CacheOfGroundTruthModelsWithCombinedTimes.Models;
        end % get.GroundTruthModelsWithCombinedTimes

        function models = get.GroundTruthModelsWithIndividualTimes(obj)
            obj.CacheOfGroundTruthModelsWithIndividualTimes = ...
                obj.multiplyModel(obj.GroundTruthModel, ...
                    obj.Configurations, ...
                    false, ... % DON'T combine times
                    false, ... % Ground-truth models, so DON'T refit
                    obj.CacheOfGroundTruthModelsWithIndividualTimes);
            models = ...
                obj.CacheOfGroundTruthModelsWithIndividualTimes.Models;           
        end % get.GroundTruthModelsWithIndividualTimes

        function models = get.GuessedModelsWithCombinedTimes(obj)
            obj.CacheOfGuessedModelsWithCombinedTimes = ...
                obj.multiplyModel(obj.GuessedModel, ...
                    obj.Configurations, ...
                    true, ... % Combine times
                    true, ... % Guessed models, so DO refit to new data
                    obj.CacheOfGuessedModelsWithCombinedTimes);
            models = ...
                obj.CacheOfGuessedModelsWithCombinedTimes.Models;
        end % get.GuessedModelsWithCombinedTimes

        function models = get.GuessedModelsWithIndividualTimes(obj)
            obj.CacheOfGuessedModelsWithIndividualTimes = ...
                obj.multiplyModel(obj.GuessedModel, ...
                    obj.Configurations, ...
                    false, ... % DON'T combine times
                    true, ... % Guessed models, so DO refit to new data
                    obj.CacheOfGuessedModelsWithIndividualTimes);
            models = ...
                obj.CacheOfGuessedModelsWithIndividualTimes.Models;
        end % get.GuessedModelsWithIndividualTimes

        function obj = MBExperimentDesigner(...
                configurations, ...
                guessedModel, ...
                initialStrategy, ...
                initialDesign, ...
                groundTruthModel, ...
                rngSeed)
            %MBExperimentDesigner Construct an instance of the experiment
            %designer class. This constructor only takes parameters for
            %those properties that must be set and remain fixed at 
            %construction time.
            %
            %In order for simulation to be possible, a ground-truth model
            %must have been specified, and since it does not make sense for
            %ground truth (from the perspective of a sequence of experiment
            %designs) to change in the course of that sequence, we require
            %that it be set at construction time if at all and prevent
            %modifications to it thereafter.
            %
            %In order for model-based experiment design to occur, a guessed
            %model must of course be provided. While that model can change
            %and likely will change over time, we also require that the
            %initial version be provided at construction time for clarity
            %of design.

            arguments
                configurations (1, :)
                guessedModel (1, 1) MBExperimentModel
                initialStrategy (1, 1) ...
                    MBAbstractExperimentDesignStrategy ...
                    {mustBeScalarOrEmpty} = []
                initialDesign MBExperimentDesign {mustBeScalarOrEmpty} = []
                groundTruthModel MBExperimentModel ...
                    {mustBeScalarOrEmpty} = []
                % RandStream bounds integer (non-shuffle) seeds in uint32
                rngSeed uint32 {mustBeScalarOrEmpty} = []
            end

            obj.Configurations = configurations;

            obj.GuessedModel = guessedModel;

            if ~isempty(initialStrategy)
                % If no initial strategy is provided, we will use the
                % default (random) strategy.

                obj.Strategy = initialStrategy;
            end

            if isempty(initialDesign)
                initialDesign = MBExperimentDesign(configurations);
                obj.NextExperimentDesign = ...
                    obj.Strategy.apportionObservations(initialDesign);
            else
                obj.NextExperimentDesign = initialDesign;
            end

            obj.GroundTruthModel = groundTruthModel;
            if isempty(obj.GroundTruthModel)
                obj.UseEmpiricalData = true;
            end

            if ~isempty(rngSeed)
                obj.LocalRNGStream = ...
                    RandStream.create("mt19937ar", "Seed", rngSeed);
            end

            obj.ReadyToPerformNextRound = true;
            % 
            % properties
            %     
            %     CumulativeNumbersOfObservations (1, :) uint64       
        end % MBExperimentDesigner

        function performNextRound(obj, pathsToEmpiricalData)
            %performNextRound performs a single round of experimentation
            %according to the most recently specified experiment design,
            %i.e., an apportionment of available observations across
            %candidate experiment configurations.
            %
            %The way the experiment round is performed depends upon the
            %mode: empirical versus simulated.
            %   Empirical: The user may specify a list of paths (absolute
            %   or relative) to files containing data to be used. All such
            %   files, IF NOT PREVIOUSLY SUPPLIED TO THE DESIGNER, will be
            %   read and appended to a table of empirical data. The subset
            %   of observations in that table NOT used by this designer in
            %   any previous round will be sampled according to the design.
            %   The guessed model will be refit according to all empirical
            %   data seen thus far, including in the nascent round.
            %
            %   Simulated: The user should not specify any paths and will
            %   be warned when having done so. A ground-truth model must
            %   have been specified and will be solved using the FSP. The
            %   solution will then be downsampled according to the design.
            %   All such data across all rounds of SIMULATED
            %   experimentation will be stored in a dedicated table.
            %   The guessed model will be refit according to all simulated
            %   data seen thus far, including in the nascent round.

            arguments
                obj (1, 1) MBExperimentDesigner
                pathsToEmpiricalData (1, :) string = []
            end

            % First, ensure that there is a next round to perform, i.e.,
            % that a new, unperformed round has been designed.

            if ~obj.ReadyToPerformNextRound
                error("performNextRound was called out of order. " + ...
                    "This means that either a suitable model and " + ...
                    "configurations have not been specified, or " + ...
                    "designNextRound has not been called since " + ...
                    "the last call to performNextRound.")
            else
                obj.PerformingRoundNumber = obj.PerformingRoundNumber + 1;
            end

            % Next, switch to using our local RNG stream, if available.
            % This is the last point before randomness will be utilized,
            % e.g., in sampling from available data. It is also the last
            % point where we are guaranteed not to error.

            if ~isempty(obj.LocalRNGStream)
                obj.GlobalRNGStream = ...
                    RandStream.setGlobalStream(obj.LocalRNGStream);
            else
                obj.GlobalRNGStream = RandStream.getGlobalStream();
            end

            try
                % Next, handle any provided empirical data files.

                if ~isempty(pathsToEmpiricalData)
                    if obj.UseEmpiricalData
                        obj.buildDataTable(pathsToEmpiricalData);
                    else % Empirical mode with data paths provided                
                        warning("Paths to empirical data files " + ...
                            "shouldn't be provided when designing " + ...
                            "experiments in simulation mode!");
                    end % Simulation mode with data paths provided
                end % Data paths provided

                % Next, simulate the system, if appropriate, according to 
                % the current design.

                if ~obj.UseEmpiricalData
                    outputFiles = obj.simulateSystem();
                    obj.buildDataTable(outputFiles);
                end % Simulation mode

                % Finally, regardless of the experimental mode, sample
                % according to the current design. This design will 
                % therefore have been used, so we are not ready to perform
                % a new round:

                obj.sampleObservationsPerDesign();

                obj.ReadyToPerformNextRound = false;

                RandStream.setGlobalStream(obj.GlobalRNGStream);
            catch ME
                % Reset the global stream and rethrow the error:
                RandStream.setGlobalStream(obj.GlobalRNGStream);
                rethrow ME
            end % [Performance, with potential faults]
        end % performNextRound        

        function set.Configurations(obj, configs)
            % This setter allows the user the flexibility to set the
            % configurations by providing either a list of configurations
            % proper (i.e., MBExperimentConfiguration) or a list of
            % CONFIGURABLES (i.e., MBAbstractExperimentConfigurable). In
            % the latter case, the configurables will be "multiplied" into
            % distinct configurations, each with a single value for each
            % configurable.

            arguments
                obj (1, 1) MBExperimentDesigner
                configs (1, :)
            end

            configurations = MBExperimentDesigner.resolveConfigurations(...
                configs);
            configurationsMustBeUnique(configurations);
            obj.Configurations = configurations;
        end % set.Configurations
    end % Public methods
end % MBExperimentDesigner

function configurationsMustBeUnique(configurations)
    arguments
        configurations (1, :) MBExperimentConfiguration
    end

    % There is no elegant way to validate uniqueness in this case other
    % than by pairwise comparison; we will exploit symmetry, however.

    for configIdx1 = 1:length(configurations)
        for configIdx2 = (1 + configIdx1):length(configurations)
            if configurations(configIdx1) == configurations(configIdx2)
                error("In a list of configurations that should " + ...
                    "have been unique, two were identical, at " + ...
                    "indices " + configIdx1 + " and " + configIdx2);
            end
        end
    end
end