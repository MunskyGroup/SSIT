classdef SSITMultiModel
    %SSITMultiModel allows for the combination of multiple SSIT models to
    %   fit a common set of parameters to multiple different data sets.

    properties
        parameterConstraints   % Function that applies prior constraints on parameters
        parameterIndices       % Set of vectors that map full parameter set to corresponding model
        SSITModels             % Set of SSIT models and their associated data sets
        fspStateSpaces         % Set of FSP StateSpaces
    end
     
    properties (Dependent)
        logLikelihoodFunctions % Set of functions that compute the log likelihood of data given each model and provided parameters
    end

    methods
        function SMM = SSITMultiModel(SSITMods,parIndices,parConstraints,stateSpaces)
            % SSITMultiModel Construct an instance of this class.
            arguments
                SSITMods = [];
                parIndices = [];
                parConstraints = [];
                stateSpaces = [];
            end   
            SMM.SSITModels=SSITMods;
            SMM.parameterIndices = parIndices;
            SMM.parameterConstraints = parConstraints;
            if ~isempty(stateSpaces)
                SMM.fspStateSpaces = stateSpaces;
            else
                nMod = length(SSITMods);
                SMM.fspStateSpaces = cell(1,nMod);
            end
        end

        function SMM = addModel(SMM,SSITMods,parIndices,parConstraints,stateSpaces)
            % Add a model or change the definition of the prior
            % constraints.
            arguments
                SMM
                SSITMods;
                parIndices;
                parConstraints = [];
                stateSpaces = [];
            end            
            for i = 1:length(SSITMods)
                SMM.SSITModels(end+1)=SSITMods(i);
                SMM.parameterIndices(end+1)=parIndices(i);
            end
            if ~isempty(stateSpaces)
                SMM.fspStateSpaces(end+1)=stateSpaces(i);
            end

            if ~isempty(parConstraints)
                SMM.parameterConstraints = parConstraints;
            end
        end

        function SMM = initializeStateSpaces(SMM,boundGuesses)
            arguments
                SMM
                boundGuesses = [];
            end
            nMod = length(SMM.SSITModels);
            for i = 1:nMod
                %% Solve the model using the FSP
                Model = SMM.SSITModels{i};
                Model.fspOptions.fspTol = 1e-4;
                if ~isempty(boundGuesses)
                    Model.fspOptions.bounds = boundGuesses{i};
                else 
                    Model.fspOptions.bounds = [];
                end
                [fspSoln,SMM.SSITModels{i}.fspOptions.bounds] = Model.solve;
                SMM.fspStateSpaces{i} = fspSoln.stateSpace;
            end
        end

        function SMM = updateModels(SMM,parameters,makeplot)
            % Updates parameters of the models to provided values and makes
            % plots of the results.
            arguments
                SMM
                parameters
                makeplot = true
            end
            Nmods = length(SMM.SSITModels);
            for i = 1:Nmods
                SMM.SSITModels{i}.parameters(:,2) = num2cell(parameters(SMM.parameterIndices{i}));
                if makeplot
                    SMM.SSITModels{i}.makeFitPlot([],5,100*i+[1:4]);
                end
            end
        end

        function objFuns = get.logLikelihoodFunctions(SMM)
            % objectiveFunction creation.  
            %     Loops through existing models and creates an objective
            %     function for each that returns the log-likelihood of the
            %     data given that model.
            Nmods = length(SMM.SSITModels);
            for i = Nmods:-1:1
                objFuns{i} = @(x)SMM.SSITModels{i}.computeLikelihood(...
                    x(SMM.parameterIndices{i}),...
                    SMM.fspStateSpaces{i});
            end
        end

        function totalLogLikelihood = computeTotalLogLikelihood(SMM,parameterGuess)
            % Method that computes the total log likeihood of all data and
            % all models for the provided parameter combination.
            Nmods = length(SMM.logLikelihoodFunctions);
            logLs = zeros(1,Nmods);
            for i = 1:Nmods
                logLs(i) = SMM.logLikelihoodFunctions{i}(parameterGuess);
            end

            totalLogLikelihood = sum(logLs);

            % Apply prior constraints.
            if ~isempty(SMM.parameterConstraints)
                totalLogLikelihood = totalLogLikelihood + SMM.parameterConstraints(parameterGuess);
            end
        end

        function [pars,likelihood] = maximizeLikelihood(SMM,parGuess,fitOptions,fitAlgorithm)
            % Search parameter space to determine which sets maximize the
            % likelihood function.  
            arguments
                SMM
                parGuess
                fitOptions = optimset('Display','iter','MaxIter',10)
                fitAlgorithm = 'fminsearch'
            end

            x0 = log(parGuess);

            objFun = @(x)-SMM.computeTotalLogLikelihood(exp(x));  % We want to MAXIMIZE the likelihood.

            switch fitAlgorithm
                case 'fminsearch'
                    [x0,likelihood]  = fminsearch(objFun,x0,fitOptions);

                case 'particleSwarm'
                    SMM.fspOptions.fspTol=inf;
                    rng('shuffle')
                    OBJps = @(x)objFun(x');
                    LB = -5*ones(size(x0'));
                    UB = 5*ones(size(x0'));
                    initSwarm = repmat(x0',fitOptions.SwarmSize-1,1);
                    initSwarm = [x0';initSwarm.*(1+0.1*randn(size(initSwarm)))];
                    fitOptions.InitialSwarmMatrix = initSwarm;
                    [x0,likelihood] = particleswarm(OBJps,length(x0),LB,UB,fitOptions);

                case 'MetropolisHastings'

                    allFitOptions.isPropDistSymmetric=true;
                    allFitOptions.thin=1;
                    allFitOptions.numberOfSamples=1000;
                    allFitOptions.burnIn=100;
                    allFitOptions.progress=true;
                    allFitOptions.proposalDistribution=@(x)x+0.01*randn(size(x));
                    allFitOptions.numChains = 1;
%                     allFitOptions.useFIMforMH = false;
%                     allFitOptions.CovFIMscale = 0.6;
                    allFitOptions.logForm = true;
                    
%                     j=1;
%                     while exist(['TMPmh_',num2str(j),'.mat'],'file')
%                         j=j+1;
%                     end                    
%                     allFitOptions.saveFile = ['TMPmh_',num2str(j),'.mat'];
                    fNames = fieldnames(fitOptions);
                    for i=1:length(fNames)
                        allFitOptions.(fNames{i}) = fitOptions.(fNames{i});
                    end

                    if allFitOptions.logForm
                        OBJmh = @(x)SMM.computeTotalLogLikelihood(exp(x));  % We want to MAXIMIZE the likelihood.
                        x0 = log(parGuess);
                    else
                        OBJmh = @(x)SMM.computeTotalLogLikelihood(x);  % We want to MAXIMIZE the likelihood.
                        x0 = (parGuess);
                    end

% Need New Method for FIM with MultiModel
%                   if allFitOptions.useFIMforMetHast
%                         TMP = SMM;
%                         TMP.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
%                         [sensSoln] = TMP.solve;  % Solve the sensitivity problem
% 
%                         fimResults = TMP.computeFIM(sensSoln.sens);
%                         [FIM] = TMP.evaluateExperiment(fimResults,TMP.dataSet.nCells);
% 
%                         if allFitOptions.logForm
%                             FIMlog = diag([TMP.parameters{SMM.fittingOptions.modelVarsToFit,2}])*...
%                                 FIM(SMM.fittingOptions.modelVarsToFit,SMM.fittingOptions.modelVarsToFit)*...
%                                 diag([TMP.parameters{SMM.fittingOptions.modelVarsToFit,2}]);
%                         else
%                             FIMlog = FIM(SMM.fittingOptions.modelVarsToFit,SMM.fittingOptions.modelVarsToFit);
%                         end
%                         
%                         covLog = FIMlog^-1;
%                         covLog = allFitOptions.CovFIMscale*(covLog+covLog')/2;
%                         allFitOptions.proposalDistribution=@(x)mvnrnd(x,covLog);                       
%                     end
                    
                    rng('shuffle')
                    if allFitOptions.numChains==1
                        [otherResults.mhSamples,otherResults.mhAcceptance,otherResults.mhValue,x0,likelihood] = ...
                            ssit.parest.metropolisHastingsSample(x0',allFitOptions.numberOfSamples,...
                            'logpdf',OBJmh,'proprnd',allFitOptions.proposalDistribution,...
                            'symmetric',allFitOptions.isPropDistSymmetric,...
                            'thin',allFitOptions.thin,'nchain',1,'burnin',allFitOptions.burnIn,...
                            'progress',allFitOptions.progress,...
                            'saveFileName',allFitOptions.saveFile);
                    else
                        try
                            parpool
                        catch
                        end
                        allFitOptions.progress=0;
                        clear tmpMH*
                        parfor iChain = 1:allFitOptions.numChains
                            [mhSamples, mhAcceptance, mhValue,xbest,fbest] = ...
                                ssit.parest.metropolisHastingsSample(x0',allFitOptions.numberOfSamples,...
                                'logpdf',OBJmh,'proprnd',allFitOptions.proposalDistribution,'symmetric',...
                                allFitOptions.isPropDistSymmetric,...
                                'thin',allFitOptions.thin,'nchain',1,'burnin',allFitOptions.burnIn,...
                                'progress',allFitOptions.progress);
                            tmpMHSamp(iChain) = {mhSamples};
                            tmpMHAcceptance(iChain) = {mhAcceptance};
                            tmpMHValue(iChain) = {mhValue};
                            tmpMHxbest(iChain) = {xbest};
                            tmpMHfbest(iChain) = fbest;
                        end
                        [~,jBest] = max(tmpMHfbest);
                        x0 = tmpMHxbest{jBest}';
                        otherResults.mhSamples = tmpMHSamp;
                        otherResults.mhAcceptance = tmpMHAcceptance;
                        otherResults.mhValue = tmpMHValue;
                        clear tmpMH*
                    end
                    % If fit was in linear space, need to convert to log
                    % space before returning parameters.
                    if ~allFitOptions.logForm
                        pars = log(x0);
                    end

            end
            pars = exp(x0);
        end
    end
end