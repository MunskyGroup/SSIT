classdef MetropolisHastingsRunner
    properties
        DataType (1, 1) ExperimentalDataType = ExperimentalDataType.Simulated
        FIMCovarianceScale (1, 1) double = 1
        FIMScale (1, 1) string {mustBeNonempty} = 'log'
        FIMToUse double = []
        IsTuning (1, 1) logical = true
        Model (1, 1) DiscoverableModel
        ModelFSPStateSpace = []
        ModelHasData (1, 1) logical = false
        NumberForBurnMH (1, 1) double {mustBeInteger, mustBePositive} = 1
        NumberForThinMH (1, 1) double {mustBeInteger, mustBePositive} = 2
        NumberOfObservations (1, 1) double {mustBeInteger, mustBePositive} = 1
        NumberOfSamplesForProduction (1, 1) double {mustBeInteger, mustBePositive} = 5000
        NumberOfSamplesForTuning (1, 1) double {mustBeInteger, mustBePositive} = 100
        Progress (1, 1) logical = false
        SaveFile (1, 1) string {mustBeNonempty} = ['TMP_MH_', randi(1000), '.mat']
        TestingMode (1, 1) logical = false
    end % Publicly mutable properties    
    
    properties (Access = private)
        ParameterFittingObjectiveMH % Function handle
        ProbabilityMH % Function handle
        ProposalDistribution % Function handle
        x0
    end % Private properties
    
    properties (SetAccess = private)
        NewParameters
        ResultsAcceptance
        ResultsSamples
        ResultsValue
    end % Privately mutable, publicly accessible properties

    methods
        function obj = run(obj)
            while obj.IsTuning
                disp(['Starting new MH chain for tuning: ', obj.SaveFile])
                delete(obj.SaveFile); % Delete old save file if it exists

                % In our implementation of the Metropolis-Hastings
                % algorithm, we pass a function to calculate the
                % probability directly rather than to calculate an "energy"
                % of the "system," from which the probability is normally
                % calculated, according to the Maxwell-Boltzmann
                % distribution. This probability is therefore the
                % likelihood (note the positive sign) of the data given the
                % parameters. On the other hand, when we seek to fit new
                % parameters to the data using the fminsearch function, we
                % provide a function that calculates the negative
                % likelihood, so that its minimization will yield the
                % maximum likelihood.
                
                obj.ProbabilityMH = @(x) obj.Model.getLikelihood(...
                    x, obj.ModelHasData, obj.ModelFSPStateSpace);

                switch obj.DataType
                    case "Simulated"
                        % Use the true (provided) FIM if we're using
                        % simulated data

                        FIM = obj.Model.evaluateExperiment(...
                            obj.FIMToUse, obj.NumberOfObservations, ...
                            obj.Model.fittingOptions.logPriorCovariance);
                    case "Empirical"
                        obj.Model.fspOptions.fspTol = 1e-8;
                        FIMResults = obj.Model.computeFIM([], obj.FIMScale);
                        FIM = obj.Model.evaluateExperiment(...
                            FIMResults, obj.NumberOfObservations, ...
                            obj.Model.fittingOptions.logPriorCovariance);
                end
                
                FIMFree = FIM(obj.Model.fittingOptions.modelVarsToFit, ...
                    obj.Model.fittingOptions.modelVarsToFit);

                if min(eig(FIMFree)) < 0.1
                    disp(...
                        ['Warning -- FIM has one or more small ' ...
                        'eigenvalues. Reducing proposal in those ' ...
                        'directions. MH convergence may be slow.'])
                    FIMFree = FIMFree + 0.1 * eye(length(FIMFree));
                end

                CovFree = FIMFree ^ -1;
                CovFree = (CovFree + CovFree') / 2;

                obj.ProposalDistribution = ...
                    @(x) mvnrnd(x, CovFree * obj.FIMCovarianceScale);

                if obj.TestingMode
                    obj.Progress = true;
                end

                obj.runMHSample();
                delete(obj.SaveFile);

                % While we are still tuning, we will not only update the
                % model if the parameters changed but also perform another
                % search starting at the new parameter set. Also, if the
                % parameters changed (significantly) during this tuning
                % iteration, or if too low an acceptance rate was 
                % encountered, we will continue tuning.

                obj.IsTuning = obj.updateModelIfParametersChanged(true);
                if ~obj.IsTuning
                    obj.FIMCovarianceScale = obj.FIMCovarianceScale / 2;
                    obj.IsTuning = true;
                end
            end % while tuning

            disp('Starting production MH chain')

            % Run MH with tuned proposal
            obj.runMHSample();
            obj.updateModelIfParametersChanged();
        end % run

        function obj = runMHSample(obj)
            if obj.IsTuning
                numberOfSamples = obj.NumberOfSamplesForTuning;
            else
                numberOfSamples = obj.NumberOfSamplesForProduction;
            end

            [obj.ResultsSamples, obj.ResultsAcceptance, ...
                obj.ResultsValue, obj.x0] = ...
                ssit.parest.metropolisHastingsSample(...
                    log(obj.NewParameters), numberOfSamples, ...
                    'logpdf', obj.ProbabilityMH, ...
                    'proprnd', obj.ProposalDistribution, ...
                    'symmetric', true, ...
                    'thin', obj.NumberForThinMH, ...
                    'nchain', 1, ...
                    'burnin', obj.NumberForBurnMH, ...
                    'progress', obj.Progress, ...
                    'saveFileName', obj.SaveFile ...                
                );
        end % runMHSample



        function obj = set.ResultsAcceptance(obj, val)
            obj.ResultsAcceptance = val;
        end
        
        function obj = set.ResultsSamples(obj, val)
            obj.ResultsSamples = val;
        end
        
        function obj = set.ResultsValue(obj, val)
            obj.ResultsValue = val;
        end
    end % Methods

    methods (Access = private)
        function [obj, parametersChanged] = updateModelIfParametersChanged(...
                obj, searchAgainIfParametersChanged)
            arguments
                obj
                searchAgainIfParametersChanged (1, 1) logical = false
            end

            parametersChanged = false;
            if max(abs(obj.NewParameters - exp(obj.x0)) ./ obj.NewParameters) > 0
                parametersChanged = true;

                % Convert best parameters to linear space
                obj.NewParameters = exp(obj.x0); 

                % Update model and its state space
                obj.Model.fspOptions.fspTol = 1e-4;
                obj.Model.parameters(obj.Model.fittingOptions.modelVarsToFit, 2) = num2cell(obj.NewParameters);
                [FSPSolution, obj.Model.fspOptions.bounds] = obj.Model.solve;
                obj.ModelFSPStateSpace = FSPSolution;

                if searchAgainIfParametersChanged
                    % Run another fminsearch starting at the new 
                    % parameter set. Redefine the objective function to
                    % calculate the new total likelihood function.

                    obj.ParameterFittingObjectiveMH = @(x) ...
                        -obj.Model.getLikelihood(...
                            x, obj.ModelHasData, obj.ModelFSPStateSpace);
                    fitOptions = optimset('Display', 'iter', ...
                        'MaxIter', 100);

                    % TODO: Pass in fitOptions from designer, so that we
                    % can get a correct value for maxFitIter.

                    obj.NewParameters = exp(fminsearch(...
                        obj.ParameterFittingObjectiveMH, ...
                        log(obj.NewParameters), fitOptions));
                    obj.Model.parameters(obj.Model.FitParameters, 2) = ...
                        num2cell(obj.NewParameters);
                    [FSPSolution, obj.Model.fspOptions.bounds] = obj.Model.solve;
                    obj.ModelFSPStateSpace = FSPSolution;
                end % if searchAgainIfParametersChanged
            end % if parametersChanged
        end % updateModelIfParametersChanged
    end % Private methods
end