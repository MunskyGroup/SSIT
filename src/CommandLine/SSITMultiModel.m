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
        end

        %% Solve and Convolve
        function [SMM1,SMM2,conv2solnTensor] = solveandconvolve(SMM1,SMM2,boundGuesses)
            arguments
                SMM1
                SMM2
                boundGuesses = [];
            end
            nMod = length(SMM1.SSITModels);
            fspSolnsSMM1 = struct(); % Allocate structure to store solutions
            fspSolnsSMM2 = struct(); % Allocate structure to store solutions
            for i = 1:nMod
                %% Solve the model using the FSP
                Model1 = SMM1.SSITModels{i};
                Model2 = SMM2.SSITModels{i};
                Model1.fspOptions.fspTol = 1e-4;
                Model2.fspOptions.fspTol = 1e-4;
                if ~isempty(boundGuesses)
                    Model1.fspOptions.bounds = boundGuesses{i};
                    Model2.fspOptions.bounds = boundGuesses{i};
                else 
                    Model1.fspOptions.bounds = [];
                    Model2.fspOptions.bounds = [];
                end

                if strcmp(Model1.solutionScheme,'FSP')
                    [fspSoln1,SMM1.SSITModels{i}.fspOptions.bounds] = Model1.solve;
                    [fspSoln2,SMM2.SSITModels{i}.fspOptions.bounds] = Model2.solve;
                    % Initialize the structure for the current model
                    fspSolnsSMM1(i).fsp = cell(numel(fspSoln1.fsp), 1); % Cell array for FSP solutions
                    fspSolnsSMM2(i).fsp = cell(numel(fspSoln2.fsp), 1); % Cell array for FSP solutions
                    for f=1:numel(fspSoln1.fsp)
                        fspSolnsSMM1(i).fsp{f} = fspSoln1.fsp{f}; 
                        fspSolnsSMM2(i).fsp{f} = fspSoln2.fsp{f};
                    end
                    SMM1.fspStateSpaces{i} = fspSoln1.stateSpace;
                    SMM2.fspStateSpaces{i} = fspSoln2.stateSpace;
                    fspSolnsSMM1(i).stateSpace = fspSoln1.stateSpace; % Store state space
                    fspSolnsSMM2(i).stateSpace = fspSoln2.stateSpace; % Store state space
                end

                %% Convolution
                for f=1:(max(numel(fspSoln1.fsp),numel(fspSoln2.fsp)))
                    % f is time point, so solution tensors are FSP probabilities across states for each time point
                    conv2solnTensor{f} = conv2(double(fspSoln1.fsp{f}.p.data),double(fspSoln2.fsp{f}.p.data));
                    figure(f)
                    contourf(log10(conv2solnTensor{f}))
                    hold on
                end
                %% check
                for g=1:f % f is number of time points 
                    fspsoln1_sptensor{g} = double(fspSoln1.fsp{g}.p.data);
                    fspsoln2_sptensor{g} = double(fspSoln2.fsp{g}.p.data);
                    figure(g)
                    subplot(1,3,1)
                    contourf(log10(fspsoln1_sptensor{g}))
                    subplot(1,3,2)
                    contourf(log10(fspsoln2_sptensor{g}))
                    subplot(1,3,3)
                    contourf(log10(conv2solnTensor{g}))
                    hold on
                end
            end
        end


        function SMM = updateModels(SMM,parameters,makeplot, fignums)
            % Updates parameters of the models to provided values and makes
            % plots of the results.
            arguments
                SMM
                parameters = []
                makeplot = false
                fignums = [];
            end

            if isempty(parameters)
                parameters = SMM.parameters;
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
            arguments
                SMM
                parameterGuess = [];
            end
            % Method that computes the total log likeihood of all data and
            % all models for the provided parameter combination.
            Nmods = length(SMM.logLikelihoodFunctions);
            logLs = zeros(1,Nmods);

            if isempty(parameterGuess)
                parameterGuess = SMM.parameters;
            end

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

        function [pars,likelihood,otherResults,SMM] = maximizeLikelihood(SMM,parGuess,fitOptions,fitAlgorithm)
            % Search parameter space to determine which sets maximize the
            % likelihood function.  
            arguments
                SMM
                parGuess = [];
                fitOptions = optimset('Display','iter','MaxIter',1000)
                fitAlgorithm = 'fminsearch'
            end

            otherResults = [];
            if isempty(parGuess)
                parGuess = SMM.parameters;
            end

            x0 = log(parGuess);

            if isfield(fitOptions,'suppressExpansion')&&fitOptions.suppressExpansion==true
                for iModel = 1:length(SMM.SSITModels)
                    oldFspTols(iModel) = SMM.SSITModels{iModel}.fspOptions.fspTol;
                    SMM.SSITModels{iModel}.fspOptions.fspTol = inf;
                end
            end

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
                            'progress',allFitOptions.progress);
                            %'saveFileName',allFitOptions.saveFile);
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
            SMM.parameters = pars;
            SMM = SMM.updateModels(pars);

            if isfield(fitOptions,'suppressExpansion')&&fitOptions.suppressExpansion==true
                for iModel = 1:length(SMM.SSITModels)
                    SMM.SSITModels{iModel}.fspOptions.fspTol = oldFspTols(iModel);
                end
            end

        end
        function compareParameters(SMM,fignum,relative)
            % This function makes a heatmap plot to compare parameters in a
            % multi model.  This only make sense when all of the sub-models
            % are the same.
            arguments
                SMM
                fignum = [];
                relative = false;
            end
            if isempty(fignum)
                figure;
            else
                figure(fignum);
            end

            nMod = length(SMM.SSITModels);
            nPars = size(SMM.SSITModels{1}.parameters,1);
            plotMat = zeros(nMod,nPars);
            for iMod = 1:nMod
                plotMat(iMod,:) = SMM.parameters(SMM.parameterIndices{iMod});
            end
            if relative
                plotMatRel = plotMat./repmat(mean(plotMat),nMod,1);
                plotMatRel(iMod+1,nPars+1) = 0;
                pcolor(log2(plotMatRel));
            else
                plotMat(iMod+1,nPars+1) = 0;
                pcolor(log10(plotMat));
            end
            set(gca,'xtick',[0.5:1:nPars+0.5],'XTickLabel',SMM.SSITModels{1}.parameters(:,1))
            set(gca,'ytick',[0.5:1:nMod+0.5],'YTickLabel',1:nMod,'FontSize',15)
            xlabel('Parameter Name')
            ylabel('Model Number')

            cb = colorbar;
            if relative
                cb.Label.String = 'log2(deviation from mean)';
            else
                cb.Label.String = 'log10(value)';
            end
        end
    end
    methods (Static)
        function SMM = createCrossValMultiModel(Model,DataFileName, ...
                LinkedSpecies,ConditionsGlobal,ConditionsReplicas, ...
                Log10Constraints, stateSpace)             
            arguments
                Model           % SSIT Model
                DataFileName    % String - name of data set
                LinkedSpecies   % String identifying which species are which columns of data file
                ConditionsGlobal% String - conserved settings for data selection
                ConditionsReplicas % String - additional settings for each replica
                Log10Constraints % Double - expected log10 uncertainty in parameters between replicas
                %                            (0=no variation, inf= no relationship)
                stateSpace =[];
            end

            nReps = length(ConditionsReplicas);

            % Make sure that modelVarsToFit is index numbers for free
            % parameters.
            if ischar(Model.fittingOptions.modelVarsToFit)&&strcmp(Model.fittingOptions.modelVarsToFit,'all')
                Model.fittingOptions.modelVarsToFit = find(ones(1,size(Model.parameters,1),'logical'));
            elseif islogical(Model.fittingOptions.modelVarsToFit)
                Model.fittingOptions.modelVarsToFit = find(Model.fittingOptions.modelVarsToFit);
            end

            nPars = length(Log10Constraints);
            if nPars~=length(Model.fittingOptions.modelVarsToFit)
                error('Length of constraints must match number of free parameters in template model.')
            end

            if ~isempty(ConditionsGlobal{1})||~isempty(ConditionsGlobal{2})
                error('Cross Validation replica creation only supported using complex conditions approach (third entry of conditions entry)')
            end

            % Match parameters in CV models to originals and to collected
            % set of all free parameters.
            indsFree = find(Log10Constraints~=0);
            parIndices = cell(1,nReps);
            parIndices{1} = (1:nPars);
            mPars = nPars;
            matchedInds = zeros(nReps,length(indsFree));
            matchedInds(1,:) = indsFree;
            for iRep = 2:nReps
                parIndices{iRep} = 1:nPars;
                parIndices{iRep}(indsFree) = mPars+1:mPars+length(indsFree);
                matchedInds(iRep,:) = parIndices{iRep}(indsFree);
                mPars = mPars+length(indsFree);
            end

            % Define penalty constraint on matched parameters.
            parConstraints = @(x)SSITMultiModel.computeReplicaMismatch(x,indsFree,matchedInds,Log10Constraints);

            SSITMods = cell(1,nReps);
            for iRep = 1:nReps
                conditionsLocal = {[],[],append(ConditionsGlobal{3},'&',ConditionsReplicas{iRep})};
                SSITMods{iRep} = Model.loadData(DataFileName, LinkedSpecies, conditionsLocal);
            end

            SMM = SSITMultiModel(SSITMods,parIndices,parConstraints,stateSpace);

            % Initialize statespaces.
            SMM = SMM.initializeStateSpaces;


        end
        function parConstraints = computeReplicaMismatch(pars,indsFree,matchedInds,Log10Constraints)
            % compute the likelihood of given parameter variation given
            % lognormal distribution with specified log10 deviation
        
            % matchedInds is nReps x nConstrainedPars; its entries are global
            % indices into 'pars', so we can index 'pars' directly:
            log10ParsMatched = log10(pars(matchedInds));   % size: nReps x nConstrainedPars
        
            nReps = size(matchedInds,1);
        
            % Mean across replicas (row-wise mean), subtract, and square
            deviation2 = sum(log10ParsMatched - repmat(mean(log10ParsMatched,1), nReps, 1)).^2;
        
            % Only penalize parameters with nonzero Log10Constraints (those in indsFree)
            parConstraints = -sum(deviation2 ./ (2 * Log10Constraints(indsFree).^2));
        end

    end
end
