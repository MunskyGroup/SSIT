classdef SSITMultiModel
    %SSITMultiModel allows for the combination of multiple SSIT models to
    %   fit a common set of parameters to multiple different data sets.

    properties
        parameterConstraints   % Function that applies prior constraints on parameters
        parameterIndices       % Set of vectors that map full parameter set to corresponding model
        SSITModels             % Set of SSIT models and their associated data sets
        fspStateSpaces         % Set of FSP StateSpaces
        FIM = [];              % FIM results
        parameters  =[];       % Parameter Values
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
            for i = 1:nMod
                SMM.parameters(1,SMM.parameterIndices{i}) = ...
                    [SMM.SSITModels{i}.parameters{SMM.SSITModels{i}.fittingOptions.modelVarsToFit,2}];
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

            nMod = length(SMM.SSITModels);
            for i = 1:nMod
                SMM.parameters(1,SMM.parameterIndices{i}) = ...
                    [SMM.SSITModels{i}.parameters{SMM.SSITModels{i}.fittingOptions.modelVarsToFit,2}];
            end
        end

        function SMM = initializeStateSpaces(SMM,boundGuesses)
            arguments
                SMM
                boundGuesses = [];
            end
            nMod = length(SMM.SSITModels);
            fspSolnsSMM = struct(); % Allocate structure to store solutions
            for i = 1:nMod
                %% Solve the model using the FSP
                Model = SMM.SSITModels{i};
                Model.fspOptions.fspTol = 1e-4;
                if ~isempty(boundGuesses)
                    Model.fspOptions.bounds = boundGuesses{i};
                else 
                    Model.fspOptions.bounds = [];
                end

                if strcmp(Model.solutionScheme,'FSP')
                    [fspSoln,SMM.SSITModels{i}.fspOptions.bounds] = Model.solve;
                    % Initialize the structure for the current model
                    fspSolnsSMM(i).fsp = cell(numel(fspSoln.fsp), 1); % Cell array for FSP solutions
                    for f=1:numel(fspSoln.fsp)
                        fspSolnsSMM(i).fsp{f} = fspSoln.fsp{f}; 
                    end
                    SMM.fspStateSpaces{i} = fspSoln.stateSpace;
                    fspSolnsSMM(i).stateSpace = fspSoln.stateSpace; % Store state space
                end
            end
            save('combinedGR_a_fspSolns.mat', 'fspSolnsSMM'); % Save fspSolns for accessibility
        end

        function SMM = updateModels(SMM,parameters,makeplot, fignums)
            % Updates parameters of the models to provided values and makes
            % plots of the results.
            arguments
                SMM
                parameters
                makeplot = true
                fignums = [];
            end
            Nmods = length(SMM.SSITModels);
            for i = 1:Nmods
                SMM.SSITModels{i}.parameters(SMM.SSITModels{i}.fittingOptions.modelVarsToFit,2) = ...
                    num2cell(parameters(SMM.parameterIndices{i}));
                if makeplot
                    solnType = SMM.SSITModels{i}.solutionScheme;
                    SMM.SSITModels{i}.solutionScheme = 'FSP';
                    if isempty(fignums)
                        SMM.SSITModels{i}.makeFitPlot([],1,100*i+[1:4]);
                    else
                        SMM.SSITModels{i}.makeFitPlot([],1,fignums(i,1:4));
                    end
                    SMM.SSITModels{i}.solutionScheme = solnType;
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
                switch SMM.SSITModels{i}.solutionScheme
                    case 'ode'
                        objFuns{i} = @(x)SMM.SSITModels{i}.computeLikelihoodODE(...
                            x(SMM.parameterIndices{i}));
                    otherwise
                        objFuns{i} = @(x)SMM.SSITModels{i}.computeLikelihood(...
                            x(SMM.parameterIndices{i}),...
                            SMM.fspStateSpaces{i});
                end
            end
        end

        function totalLogLikelihood = computeTotalLogLikelihood(SMM,parameterGuess)
            % Method that computes the total log likeihood of all data and
            % all models for the provided parameter combination.
            Nmods = length(SMM.logLikelihoodFunctions);
            logLs = zeros(1,Nmods);

            G = cell(1,Nmods);
            for i = 1:Nmods
                G{i} = SMM.logLikelihoodFunctions{i};
            end
            for i = 1:Nmods
                logLs(i) = G{i}(parameterGuess);
            end

            totalLogLikelihood = sum(logLs);

            % Apply prior constraints.
            if ~isempty(SMM.parameterConstraints)
                totalLogLikelihood = totalLogLikelihood + SMM.parameterConstraints(parameterGuess(:));
            end
        end

        function SMM = computeFIMs(SMM,sensSoln,scale,MHSamples)
            arguments
                SMM
                sensSoln = [];
                scale = 'lin';
                MHSamples = [];
            end
            % This function will compute the individual FIM matrices for
            % the time points in the current experiment design. 
            Nmods = length(SMM.SSITModels);
            FIMlocal.fims = cell(1,Nmods);
            
            uniquePars = [];
            for i = 1:Nmods
                uniquePars = sort(unique([uniquePars,SMM.parameterIndices{i}]));
            end            
            
            FIMlocal.totalFIM = zeros(max(uniquePars));

            for i = 1:Nmods
                if ~isempty(MHSamples)
                    MHSamplesMod =  MHSamples(:,SMM.parameterIndices{i});
                else
                    MHSamplesMod = MHSamples;
                end
                FIMlocal.fims{i} = SMM.SSITModels{i}.computeFIM(sensSoln,scale,MHSamplesMod);
               
                J = SMM.SSITModels{i}.fittingOptions.modelVarsToFit;
                for iT = 1:length(SMM.SSITModels{i}.tSpan)
                    FIMlocal.totalFIM(SMM.parameterIndices{i},SMM.parameterIndices{i}) = ...
                        FIMlocal.totalFIM(SMM.parameterIndices{i},SMM.parameterIndices{i}) +...
                        SMM.SSITModels{i}.dataSet.nCells(iT)*FIMlocal.fims{i}{iT}(J,J);
                end
            end
            SMM.FIM = FIMlocal;
        end

        function [pars,likelihood,otherResults] = maximizeLikelihood(SMM,parGuess,fitOptions,fitAlgorithm)
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
                    allFitOptions.useFIMforMetHast = false;
                    allFitOptions.CovFIMscale = 0.6;
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
                        OBJmh = @(x)SMM.computeTotalLogLikelihood(exp(x'));  % We want to MAXIMIZE the likelihood.
                        x0 = log(parGuess)';
                    else
                        OBJmh = @(x)SMM.computeTotalLogLikelihood(x');  % We want to MAXIMIZE the likelihood.
                        x0 = (parGuess)';
                    end

                    if allFitOptions.useFIMforMetHast
                        TMP = SMM;
                        if isempty(TMP.FIM)
                            TMP = TMP.computeFIMs;
                        end
                        FIMLocal = TMP.FIM.totalFIM;

                        if allFitOptions.logForm
                            FIMlog = diag(TMP.parameters)*...
                                FIMLocal*...
                                diag(TMP.parameters);
                        else
                            FIMlog = FIMLocal(SMM.fittingOptions.modelVarsToFit,SMM.fittingOptions.modelVarsToFit);
                        end

                        if allFitOptions.logForm&&min(eig(FIMlog))<1
                            disp('Warning -- FIM has one or more small eigenvalues.  Reducing proposal width to 1x in those directions. MH Convergence may be slow.')
                            FIMlog = FIMlog + 1*eye(length(FIMlog));
                        end

                        covLog = FIMlog^-1;
                        covLog = allFitOptions.CovFIMscale*(covLog+covLog')/2;
                        allFitOptions.proposalDistribution=@(x)mvnrnd(x,covLog);
                    end
                    
                    rng('shuffle')
                    if allFitOptions.numChains==1
                        [otherResults.mhSamples,otherResults.mhAcceptance,otherResults.mhValue,x0,likelihood] = ...
                            ssit.parest.metropolisHastingsSample(x0',allFitOptions.numberOfSamples,...
                            'logpdf',OBJmh,'proprnd',allFitOptions.proposalDistribution,...
                            'symmetric',allFitOptions.isPropDistSymmetric,...
                            'thin',allFitOptions.thin,'nchain',1,'burnin',allFitOptions.burnIn,...
                            'progress',allFitOptions.progress,...
                            'saveFileName',allFitOptions.saveFile);
                        % delete(allFitOptions.saveFile);
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
