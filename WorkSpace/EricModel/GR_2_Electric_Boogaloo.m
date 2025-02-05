%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations

%% TEST NEW GR MODELS

close all 
clear
addpath(genpath('../../src'));
addpath('tmpPropensityFunctions');
loadPrevious = false;
savedWorkspace = 'workspace_Boogaloo';

fitOptions = optimset('Display','iter','MaxIter',300);
 
%%
GR = input('(0) FSP for GR-alpha only (base model);\n(1) Convolve GR-alpha + GR-beta;\n(2) ODES for GR-alpha + GR-beta.\nChoose your destiny: ');

switch GR
    case 0
    % GR-alpha setup    
        disp('You have chosen the base model, GR-alpha.')
        %% STEP 0.B. -- Create Base Model for GR Only
        % Here, I set up the model for the GR translocation dynamics.            
        % Create blank SSIT model.
        ModelGR = SSIT;
    
        % Set species names.
        ModelGR.species = {'cytGR';'nucGR'};
    
        % Set initial condition.
        ModelGR.initialCondition = [20;1];
    
        % Define propensity functions and input signals:
        ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';...
            'kg1';'gg1*cytGR';'gg2*nucGR'};
    
        % Define stoichiometry.
        ModelGR.stoichiometry = [-1,1,1,-1,0;...
            1,-1,0,0,-1];
    
        % Specify parameter guesses.
        ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
            'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',14e-5;...
            'gg1',1e-5;'gg2',1e-6;'MDex',5;'Dex0',100});
    
        % Print visual summary of the model
        ModelGR.summarizeModel
    
        %% The log prior will be applied to the fit to multiple models as an additional constraint.
        log10PriorMean = [-1 -1 0 -2,...
            -1 -3 -2 -1 -2 -2 -2 0.5 2];
        log10PriorStd = 2*ones(1,13);
        
        % So it is left out of the prior, since we only want it to be calculated once.
        ModelGR.fittingOptions.logPrior = [];  
    
        
        ModelGR.fspOptions.initApproxSS = true;
        ModelGR.fittingOptions.modelVarsToFit = (5:12);
        
        ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        ModelGR = ModelGR.formPropensitiesGeneral('EricRonModGR');
        ModelGR.customConstraintFuns = {'cytGR+nucGR'};
        % TODO - Alex - Constraints for removing stiff dimensions.
        
        [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve;
        [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve(FSPGrSoln.stateSpace);
    
        % STEP 0.B.2. -- Define GR parameters
        GRpars = cell2mat(ModelGR.parameters(5:12,2))';  
    
        % STEP 0.B.3. -- Associate GR Data with Different Instances of Model (10,100nm Dex)
        GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
            '10','10',102,'GR Fit (10nM Dex)';...
            '100','100',103,'GR Fit (100nM Dex)'};
        ModelGRparameterMap = cell(1,size(GRfitCases,1));
        ModelGRfit = cell(1,size(GRfitCases,1));
        % ModelGRODEfit = cell(1,size(GRfitCases,1));
        for i=1:3
            ModelGRfit{i} = ModelGR.loadData("EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR','normgrnuc';'cytGR','normgrcyt'},...
                {'Dex_Conc',GRfitCases{i,2}});
            ModelGRfit{i}.parameters(13,:) = {'Dex0', str2num(GRfitCases{i,1})};
            ModelGRparameterMap(i) = {(1:8)};
            % parameters 1 - 8 refer to the parameter set that is relevant to
            % the entire class of models.  In this case, these refer to
            % global parameters 5:13 (GR parameters).
        end
        
        % STEP 0.B.4. -- Make Guesses for the FSP bounds
        % This is sometimes necessary when using an uninduced steady state as the
        % initial condition. You need to guess a reasonable statespace or the
        % computation of the SS can be inaccurate.
        ModelGR.customConstraintFuns = {'cytGR+nucGR'};
        for i = 1:3
            boundGuesses{i} = [0;0;30;30;30];
            % First N are lower bounds.  Next N is upper bound.  Remaining are
            % custom.
        end
    
    %% STEP 1 -- Fit GR Models.  
    % STEP 1 will need to be rerun until satisfied.  Use fitMHiters as needed.
    % TODO: Automate with statistics.
    % Set for STEP1 -- Fit GR Models
    fitIters = 3;
    fitMHiters = 2;
    
    for GR = 1:fitMHiters
        % STEP 1.A. -- Specify dataset time points.    
        for i = 1:3
            ModelGRfit{i}.tSpan = ModelGRfit{i}.dataSet.times;
        end
    
        % STEP 1.B. -- Specify log prior (NOTE: must transpose due to Matlab update that
        %     no longer correctly assumes format when adding single value vector to
        %     column vector).
    
        logPriorGR = @(x)-sum((log10(x)-log10PriorMean(5:12)').^2./(2*log10PriorStd(5:12)'.^2));
    
        % STEP 1.C. -- Combine all three GR models and fit using a single parameter set.
        for jj = 1:fitIters
            combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap,logPriorGR);
            combinedGRModel = combinedGRModel.initializeStateSpaces(boundGuesses);
            combinedGRModel = combinedGRModel.updateModels(GRpars,false);
            GRpars = combinedGRModel.maximizeLikelihood(...
                GRpars, fitOptions);
            save('EricModel_MMDex','GRpars') 
        end
    
        save('EricModelGR_MMDex','GRpars','combinedGRModel', 'ModelGRfit', 'log10PriorStd') 
    
        %% STEP 1.D. -- Compute FIM for GR parameters.
        combinedGRModel = combinedGRModel.computeFIMs([],'log');
        fimGR_withPrior = combinedGRModel.FIM.totalFIM+... % the FIM in log space.
            diag(1./(log10PriorStd(ModelGR.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.
    
        if min(eig(fimGR_withPrior))<1
            disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
            fimGR_withPrior = fimGR_withPrior + 1*eye(length(fimGR_withPrior));
        end
        covFree = fimGR_withPrior^-1;
        covFree = 0.5*(covFree+covFree');
    
        %% STEP 1.E. -- Run MH on GR Models.
        %GRpars = GRpars';
        MHFitOptions.thin=1;
        MHFitOptions.numberOfSamples=3000;
        MHFitOptions.burnIn=1000;
        MHFitOptions.progress=true;
        MHFitOptions.numChains = 1;
    
        % Use FIM computed above rather than making SSIT call 'useFIMforMetHast'
        % which forces SSIT.m to compute it within (no prior, etc.)
        MHFitOptions.useFIMforMetHast = false;
        MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
    
        MHFitOptions.saveFile = 'TMPEricMHGR.mat';
        [~,~,MHResultsGR] = combinedGRModel.maximizeLikelihood(...
            GRpars, MHFitOptions, 'MetropolisHastings');
        %delete(MHFitOptions.saveFile)
        %MHResultsGR
        %%
        figNew = figure;
        ModelGR.plotMHResults(MHResultsGR,[],'log',[],figNew)
        for i = 1:7
            for j = i+1:7
                subplot(7,7,(i-1)*7+j-1)
                CH = get(gca,'Children');
                CH(1).Color=[1,0,1]; %
                CH(1).LineWidth = 3;
            end
        end
    end 
    
    %%     STEP 1.F. -- Make Plots of GR Fit Results
    makeGRPlots(combinedGRModel,GRpars)
    
    save('EricModelGR_MMDex','GRpars','combinedGRModel','MHResultsGR') 
    save('workspace_Feb4_2024.mat','GRpars', 'ModelGRfit', 'combinedGRModel','MHResultsGR', 'log10PriorStd')

    %% STEP 2.A.1. -- Add DUSP1 to the existing GR model.
    % Copy parameters from the 100nM Dex stim case in GR.
    fitIters = 3;
    ModelGRDusp100nM = ModelGRfit{3};
    ModelGRDusp100nM.species = {'offGene';'onGene';'cytGR';'nucGR';'rna'};
    ModelGRDusp100nM.initialCondition = [2;0;24;1;5];
    ModelGRDusp100nM.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
        '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR';...
        'kr*onGene';'gr*rna'};
    ModelGRDusp100nM.stoichiometry = [-1,1,0,0,0,0,0,0,0;...
        1,-1,0,0,0,0,0,0,0;...
        0,0,-1,1,1,-1,0,0,0;...
        0,0,1,-1,0,0,-1,0,0;...
        0,0,0,0,0,0,0,1,-1];
    ModelGRDusp100nM.useHybrid = true;
    ModelGRDusp100nM.hybridOptions.upstreamODEs = {'cytGR','nucGR'};
    ModelGRDusp100nM.solutionScheme = 'FSP';
    ModelGRDusp100nM.customConstraintFuns = [];
    ModelGRDusp100nM.fspOptions.bounds = [0;0;0;2;2;400];
    ModelGRDusp100nM.fittingOptions.modelVarsToFit = 1:4;
    ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('EricModDusp1');
    duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(1:4)).^2./(2*log10PriorStd(1:4).^2));
    ModelGRDusp100nM.fittingOptions.logPrior = duspLogPrior;
    
    %% STEP 2.A.2. -- Load pre-fit parameters into model.
    if loadPrevious
        load('EricModelDusp1_MMDex','DUSP1pars')
    else
        %% Pull the DUSP1 parameters from the Model
        % Find the indices of the desired parameter names
        knc_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'knc'));
        kg1_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'kg1'));
        gg1_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'gg1'));
        gg2_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'gg2'));

        % Extract the values and store them in DUSP1pars
        DUSP1pars = [ModelGRDusp100nM.parameters{knc_idx, 2}, ...
                    ModelGRDusp100nM.parameters{kg1_idx, 2}, ...
                    ModelGRDusp100nM.parameters{gg1_idx, 2}, ...
                    ModelGRDusp100nM.parameters{gg2_idx, 2}]; 
    end
    ModelGRDusp100nM.parameters(ModelGRDusp100nM.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
    %% STEP 2.A.3. -- Load and Associate with DUSP1 smFISH Data (100nM Dex Only)
    % The commented code below would be needed to fit multiple conditions,
    % but that is not used in this case.  It is left here in case it is
    % needed in later stages of the project.

    if ~loadPrevious
        Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
        ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
        ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
        for i = 1:size(Dusp1FitCases,1)
             ModelDusp1Fit{i} = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
                 {'rna','totalNucRNA'},...
                 {'Dex_Conc','100'});
             ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
             ModelDusp1parameterMap{i} = (1:4);
             % Set Dex concentration.
             ModelDusp1Fit{i}.parameters{13,2} = str2num(Dusp1FitCases{i,1});
             ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
        end
        DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
    end

    ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','totalNucRNA'},{'Dex_Conc','100'});
    DUSP1pars = [ModelGRDusp100nM.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];

    %%    STEP 2.B. -- Fit DUSP1 model at 100nM Dex.  
    % This will likely need to be rerun after Metropolis-Hastings Search (STEP 2.D.3).
    
    % Use fitMHiters to run until satisfied. TODO: Automate by computing statistics.
    fitMHiters = 2;
    
    for DS = 1:fitMHiters
        for i = 1:fitIters
            fitOptions.suppressFSPExpansion = true;
            DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
                DUSP1pars, fitOptions);
            ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
            save('EricModelDusp1_MMDex','GRpars','DUSP1pars') 
        end
    
        %%    STEP 2.C. -- Plot predictions for other Dex concentrations.
        showCases = [1,1,1,1];
        makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1pars,Dusp1FitCases,showCases)
    
        %% STEP 2.D. -- Sample uncertainty for Dusp1 Parameters
        %   STEP 2.D.1. -- Compute sensitivity of the FSP solution
        ModelGRDusp100nM.solutionScheme = 'fspSens';
        sensSoln = ModelGRDusp100nM.solve();
        ModelGRDusp100nM.solutionScheme = 'FSP';
        %   STEP 2.D.2. -- Compute FIM
        %       define which species in model are not observed.
        ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};
        % compute the FIM
        fimResults = ModelGRDusp100nM.computeFIM(sensSoln.sens,'log');
        % In the following, the log-prior is used as a prior co-variance matrix.
        % This will be used in the FIM calculation as an FIM without new evidence 
        % being set equal to the inverse of this covariance matrix.  More rigorous
        % justification is needed to support this heuristic.
        fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
            diag(log10PriorStd(1:13).^2));
        FIMfree = fimTotal{1}(1:4,1:4);
        if min(eig(FIMfree))<1
            disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
            FIMfree = FIMfree + 1*eye(length(FIMfree));
        end
        covFree = FIMfree^-1;
        covFree = 0.5*(covFree+covFree');
    
        %%      STEP 2.D.3. -- Run Metropolis Hastings Search
        if loadPrevious
            MHDusp1File = 'MHDusp1_Dec92024';
            load(MHDusp1File)
        else
            MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
            MHFitOptions.thin=1;
            MHFitOptions.numberOfSamples=2000;
            MHFitOptions.burnIn=0;
            MHFitOptions.progress=true;
            MHFitOptions.numChains = 1;
            MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
            [DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
                [], MHFitOptions, 'MetropolisHastings');
            delete('TMPEricMHDusp1.mat')
            ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
        end
        
        save('workspaceDec9_2024.mat','ModelGRDusp100nM','DUSP1pars','fimTotal','sensSoln','combinedGRModel','MHResultsDusp1')
    
        %%      STEP 2.D.4. -- Plot the MH results
        figNew = figure;
        ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
        for j = 1:3
            for k = i:3
                subplot(3,3,(j-1)*3+k)
                CH = get(gca,'Children');
                CH(1).Color=[1,0,1]; %
                CH(1).LineWidth = 3;
            end
        end
    end
    
    %% Save results
    varNames = unique({'ModelGR'
        'GRfitCases'
        'log10PriorMean'
        'log10PriorStd'
        'GRpars'
        'MHResultsGR'
        'ModelGRparameterMap'
        'ModelGRfit'
        'boundGuesses'
        'ModelGRDusp100nM'
        'GRfitCases'
        'log10PriorMean'
        'log10PriorStd'
        'duspLogPrior'
        'DUSP1pars'
        'ModelGRfit'
        'fimResults'
        'MHResultsDusp1'
        'sensSoln'
        });
    
    save('workspace_Feb4_2024',varNames{:}) % WARNING: THIS OVERWRITE THE PREVIOUSLY SAVED WORKSPACE - TODO: FIX

    case 1
    %% GR-alpha + GR-beta setup
        disp('You have chosen GR-alpha + GR-beta by convolution.  WARNING: UNDER CONSTRUCTION')
        ModelGR_a = SSIT;
        ModelGR_b = SSIT;
        ModelGR_a.species = {'cytGR_a';'nucGR_a'};
        ModelGR_b.species = {'cytGR_b';'nucGR_b'};
        ModelGR_a.initialCondition = [5;1];
        ModelGR_b.initialCondition = [5;1];
        ModelGR_a.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'dc * cytGR_a'; 'kn2c * nucGR_a'; 'dn * nucGR_a'};
        ModelGR_b.propensityFunctions = {'kb1 * cytGR_b'; 'bb1'; 'dc * cytGR_b'; 'kn2c * nucGR_b'; 'dn * nucGR_b'};
        ModelGR_a.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        ModelGR_b.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        ModelGR_a.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'ba1',14e-5;...
                                'MDex',5;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;'Dex0',100});
        ModelGR_b.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'kb1',0.01;'bb1',14e-5});
        ModelGR_a.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        disp('The GR-alpha model is: ')
        ModelGR_a.summarizeModel
        disp('and the GR-beta model is: ')
        ModelGR_b.summarizeModel

        ModelGR_a.fspOptions.initApproxSS = true;
        ModelGR_b.fspOptions.initApproxSS = true;

        ModelGR_a.fittingOptions.modelVarsToFit = (1:8);
        ModelGR_b.fittingOptions.modelVarsToFit = (1:5);
         
        ModelGR_a = ModelGR_a.formPropensitiesGeneral('GR_a');
        ModelGR_b = ModelGR_b.formPropensitiesGeneral('GR_b');

        ModelGR_a.customConstraintFuns = {'cytGR_a+nucGR_a'};
        ModelGR_b.customConstraintFuns = {'cytGR_b+nucGR_b'};
        % % TODO - Alex - Constraints for removing stiff dimensions.

        %% Compute FSP solutions
        [FSPsoln_a,ModelGR_a.fspOptions.bounds] = ModelGR_a.solve;
        [FSPsoln_a,ModelGR_a.fspOptions.bounds] = ModelGR_a.solve(FSPsoln_a.stateSpace);

        [FSPsoln_b,ModelGR_b.fspOptions.bounds] = ModelGR_b.solve;
        [FSPsoln_b,ModelGR_b.fspOptions.bounds] = ModelGR_b.solve(FSPsoln_b.stateSpace);
        
        %% Convolution
        for f=1:(max(numel(FSPsoln_a.fsp),numel(FSPsoln_b.fsp)))
            % f is time point, so solution tensors are FSP probabilities across states for each time point
            conv2solnTensor{f} = conv2(double(FSPsoln_a.fsp{f}.p.data),double(FSPsoln_b.fsp{f}.p.data));
            figure(f)
            contourf(log10(conv2solnTensor{f}))
            hold on
        end
        %% check
        for g=1:f % f is number of time points 
            fspsoln_sptensor_a{g} = double(FSPsoln_a.fsp{g}.p.data);
            fspsoln_sptensor_b{g} = double(FSPsoln_b.fsp{g}.p.data);
            figure(g)
            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a{g}))
            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b{g}))
            subplot(1,3,3)
            contourf(log10(conv2solnTensor{g}))
            hold on
        end
         
        %% TODO: clean up scattered data on contour plots

        data = readtable('../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv');

        % Filter the data if desired for particular Time_index and/or Dex_Conc
        filteredData = data(data.time == 180 & data.Dex_Conc == 100, :);

        % Randomly sample 1000 points from the filtered data
        numSamples = min(1000, height(filteredData)); % Ensure we don’t sample more than available data
        randomIndices = randperm(height(filteredData), numSamples);
        normgrnuc = filteredData.normgrnuc(randomIndices);
        normgrcyt = filteredData.normgrcyt(randomIndices);

        % Plot contourf figures and add scatter points
        for g=1:f
            fspsoln_sptensor_a{g} = double(FSPsoln_a.fsp{g}.p.data);
            fspsoln_sptensor_b{g} = double(FSPsoln_b.fsp{g}.p.data);
    
            figure(g)
    
            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a{g}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b{g}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,3)
            contourf(log10(conv2solnTensor{g}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
        end
                

        %% Fit GR_a to 
        %% Define GR parameters
        GRpars_a = cell2mat(ModelGR_a.parameters(1:8,2))'; 
        GRpars_b = cell2mat(ModelGR_b.parameters(1:5,2))';

        % The log prior will be applied to the fit to multiple models as an additional constraint.
        log10PriorMean_a = [-2 -5 -6 -5 0.5 -3 -3 -2];
        log10PriorMean_b = [-2 -5 -6 -2 -5];
        log10PriorStd_a = 2*ones(1,8);
        log10PriorStd_b = 2*ones(1,5);
         
        % So it is left out of the prior, since we only want it to be calculated once.
        ModelGR_a.fittingOptions.logPrior = [];  
        ModelGR_b.fittingOptions.logPrior = [];          
    
        %% Associate GR Data with Different Instances of Model (1,10,100nm Dex)
        GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
            '10','10',102,'GR Fit (10nM Dex)';...
            '100','100',103,'GR Fit (100nM Dex)'};

        ModelGRparameterMap_a = cell(1,size(GRfitCases,1));
        ModelGRfit_a = cell(1,size(GRfitCases,1));
        ModelGRfit_b = cell(1,1);

        % Load data for GR-beta
        ModelGRfit_b{1} = ModelGR_b.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR_b','normgrnuc';'cytGR_b','normgrcyt'},{'Dex_Conc','100'});

        % ModelGRODEfit = cell(1,size(GRfitCases,1));
        
        % Load data, fit 3 Dex concentrations to GR-alpha's 'Dex0' parameter,
        % and set solution scheme to FSP.
        for i=1:3
            ModelGRfit_a{i} = ModelGR_a.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR_a','normgrnuc';'cytGR_a','normgrcyt'},...
                {'Dex_Conc',GRfitCases{i,2}});
            ModelGRfit_a{i}.parameters(9,:) = {'Dex0', str2num(GRfitCases{i,1})};
            ModelGRparameterMap_a(i) = {(1:8)}; % Parameters 1:8 refers to all ModelGR_a parameters with the exception of Dex0. 
        end
 
        % Make Guesses for the FSP bounds
        % This is sometimes necessary when using an uninduced steady state as the
        % initial condition. You need to guess a reasonable statespace or the
        % computation of the SS can be inaccurate.
        for i = 1:3
            boundGuesses{i} = [0;0;30;30;30];
            % First N are lower bounds.  Next N is upper bound.  Remaining are
            % custom.
        end

    % Specify log prior (NOTE: must transpose due to Matlab update that
    %     no longer correctly assumes format when adding single value vector to
    %     column vector).

    logPriorGR_a = @(x)-sum((log10(x)-log10PriorMean_a(1:8)').^2./(2*log10PriorStd_a(1:8)'.^2));
    logPriorGR_b = @(x)-sum((log10(x)-log10PriorMean_b(1:5)').^2./(2*log10PriorStd_b(1:5)'.^2));
    
    
    %% Fit GR Models.  
    % This step will need to be rerun until satisfied.  Use fitMHiters as needed.
    % TODO: Automate with statistics.

    fitIters = 3; % 3 concentrations of Dex (1,10,100nM)
    fitMHiters = 2; % Adjust as needed
    
    for GR = 1:fitMHiters
        % Specify dataset time points for GR-alpha.    
        for i = 1:3
            ModelGRfit_a{i}.tSpan = ModelGRfit_a{i}.dataSet.times;
        end

        % Specify dataset time points for GR-beta.
        ModelGRfit_b{1}.tSpan = ModelGRfit_b{1}.dataSet.times;

        %% Solve for combined GR-alpha model for three Dex_Conc=Dex0 and fit using a single parameter set.
        for jj = 1:fitIters
            % Solve
            combinedGRModel_a = SSITMultiModel(ModelGRfit_a,ModelGRparameterMap_a,logPriorGR_a);
            combinedGRModel_a = combinedGRModel_a.initializeStateSpaces(boundGuesses);
            combinedGRModel_a = combinedGRModel_a.updateModels(GRpars_a,false);
            GRpars_a = combinedGRModel_a.maximizeLikelihood(GRpars_a, fitOptions);
            save('combinedGRModel_a','GRpars_a') 
        end
        
        %% Solve for GR-beta model.
        ModelGRfit_b{1}.fspOptions.fspTol = 1e-4;
        ModelGRfit_b{1}.fspOptions.bounds = boundGuesses{1};
        [GR_b_fspSoln,ModelGRfit_b{1}.fspOptions.bounds] = ModelGRfit_b{1}.solve;
        [GR_b_fspSoln,ModelGRfit_b{1}.fspOptions.bounds] = ModelGRfit_b{1}.solve(GR_b_fspSoln.stateSpace);
        GRpars_b = ModelGRfit_b{1}.maximizeLikelihood(GRpars_b, fitOptions);
        save('ModelGRfit_b','GRpars_b','GR_b_fspSoln')             

        save('GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln')
        

        %% Compute FIM 
        % for ModelGR_a parameters.
        combinedGRModel_a = combinedGRModel_a.computeFIMs([],'log');

        fimGR_a_withPrior = combinedGRModel_a.FIM.totalFIM+... % the FIM in log space.
            diag(1./(log10PriorStd_a(ModelGR_a.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.

        % for ModelGR_b parameters.
        ModelGR_b_fimResults = ModelGRfit_b{1}.computeFIM([],'log'); % Compute individual FIMs

        fimGR_b_withPrior = ModelGRfit_b{1}.totalFim(ModelGR_b_fimResults,ModelGRfit_b{1}.dataSet.nCells,diag(1./(log10PriorStd_b(ModelGRfit_b{1}.fittingOptions.modelVarsToFit)*log(10)).^2));  % Add prior in log space.

        %
        if min(eig(fimGR_a_withPrior))<1
            disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
            fimGR_a_withPrior = fimGR_a_withPrior + 1*eye(length(fimGR_a_withPrior));
        end
        fimGR_a_covFree = fimGR_a_withPrior^-1;
        fimGR_a_covFree = 0.5*(fimGR_a_covFree+fimGR_a_covFree');

        if min(eig(fimGR_b_withPrior{1}))<1
            disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
            fimGR_b_withPrior{1} = fimGR_b_withPrior{1} + 1*eye(length(fimGR_b_withPrior{1}));
        end
        fimGR_b_covFree = fimGR_b_withPrior{1}^-1;
        fimGR_b_covFree = 0.5*(fimGR_b_covFree+fimGR_b_covFree');

    
        %% Run MH on GR Models.
        %GRpars = GRpars';
        % GR_alpha
        MHfitOptions_a.thin=1;
        MHfitOptions_a.numberOfSamples=1000;
        MHfitOptions_a.burnIn=100;
        MHfitOptions_a.progress=true;
        MHfitOptions_a.numChains = 1;
    
        % Use FIM computed above rather than making SSIT call 'useFIMforMetHast'
        % which forces SSIT.m to compute it within (no prior, etc.)
        MHfitOptions_a.useFIMforMetHast = false;
        MHfitOptions_a.proposalDistribution=@(x)mvnrnd(x,fimGR_a_covFree);
    
        MHfitOptions_a.saveFile = 'TMPEricMHGR_a.mat';
        [~,~,MHResultsGR_a] = combinedGRModel_a.maximizeLikelihood(...
            GRpars_a, MHfitOptions_a, 'MetropolisHastings');
        %delete(MHFitOptions.saveFile)
        %MHResultsGR
        %
        figNew = figure;
        ModelGR_a.plotMHResults(MHResultsGR_a,[],'log',[],figNew)
        for i = 1:7
            for j = i+1:7
                subplot(7,7,(i-1)*7+j-1)
                CH = get(gca,'Children');
                CH(1).Color=[1,0,1]; %
                CH(1).LineWidth = 3;
            end
        end
        %% GR_beta
        % MHFitOptions_b.thin=1;
        % MHFitOptions_b.numberOfSamples=1000;
        % MHFitOptions_b.burnIn=100;
        % MHFitOptions_b.progress=true;
        % MHFitOptions_b.numChains = 1;
        % 
        % % Use FIM computed above rather than making SSIT call 'useFIMforMetHast'
        % % which forces SSIT.m to compute it within (no prior, etc.)
        % MHFitOptions_b.useFIMforMetHast = false;
        % MHFitOptions_b.proposalDistribution=@(y)mvnrnd(y,fimGR_b_covFree);
        % 
        % MHFitOptions_b.saveFile = 'TMPEricMHGR_b.mat';
        % [~,~,MHResultsGR_b] = ModelGRfit_b{1}.maximizeLikelihood(...
        %     GRpars_b, MHFitOptions_b, 'MetropolisHastings');
        % %delete(MHFitOptions.saveFile)
        % %MHResultsGR
        % %
        % figNew = figure;
        % ModelGR_b.plotMHResults(MHResultsGR_b,[],'log',[],figNew)
        % for i = 1:7
        %     for j = i+1:7
        %         subplot(7,7,(i-1)*7+j-1)
        %         CH = get(gca,'Children');
        %         CH(1).Color=[1,0,1]; %
        %         CH(1).LineWidth = 3;
        %     end
        % end
    end

    %% Convolution, contour plots, and scattered data
    load('combinedGR_a_fspSolns.mat')
    data = readtable('../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv'); 

    tspan = [0,10,30,50,75,150,180];

    for model=1:size(GRfitCases,1)
        dataDex = data(data.Dex_Conc == str2double(GRfitCases{model}), :);
        for t=1:max(numel(fspSolnsSMM(model).fsp),numel(GR_b_fspSoln.fsp))
            % Figure index
            figIndex = (model - 1) * 7 + t;
            % Filter data           
            dataDexTime = dataDex(dataDex.time == tspan(t),:);

            % Randomly sample 1000 points from the filtered data
            numSamples = min(1000, height(dataDexTime)); % Ensure we don’t sample more than available data
            % check
            if numSamples==0
                GRfitCases{model}
                tspan(t)
            end
            randomIndices = randperm(height(dataDexTime), numSamples);
            normgrnuc = dataDexTime.normgrnuc(randomIndices);
            normgrcyt = dataDexTime.normgrcyt(randomIndices);
                
            % t is time point, so solution tensors are FSP probabilities across states for each time point
            conv2solnTensor_postData{t} = conv2(double(fspSolnsSMM(model).fsp{t}.p.data),double(GR_b_fspSoln.fsp{t}.p.data));
            figure(figIndex)
            contourf(log10(conv2solnTensor_postData{t}))
            hold on 
            scatter(normgrcyt, normgrnuc, 'r', 'filled')

            % Add subplots (can be commented out)
            fspsoln_sptensor_a_postData{t} = double(fspSolnsSMM(model).fsp{t}.p.data);
            fspsoln_sptensor_b_postData{t} = double(GR_b_fspSoln.fsp{t}.p.data);

            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a_postData{t}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b_postData{t}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,3)
            contourf(log10(conv2solnTensor_postData{t}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
        end
    end

    save('conv2solnTensor_postData','GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln', 'fspSolnsSMM')

    dusp1 = input('(0) GR-beta has no effect on DUSP1;\n(1) GR-beta turns off the DUSP1 gene');

    switch dusp1
        case 0
        %% STEP 2.A.1. -- Add DUSP1 to the existing GR model.
            % Copy parameters from the 100nM Dex stim case in GR.
            fitIters = 3;
            ModelGRDusp100nM = ModelGRfit_a{3};
            ModelGRDusp100nM.species = {'offGene';'onGene';'cytGR_a';'nucGR_a';
                'rna';'cytGR_b';'nucGR_b'};
            ModelGRDusp100nM.initialCondition = [2;0;24;1;5;12;12];
            ModelGRDusp100nM.propensityFunctions = {'kon*offGene*nucGR_a';'koff*onGene';
                '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex))*cytGR_a';'knc*nucGR_a';
                'kg1';'gg1*cytGR_a';'gg2*nucGR_a';'kr*onGene';'gr*rna';
                'kcn2*cyt_b';'knc*nucGR_b';'kg2';'gg1*cytGR_b';'gg2*nucGR_b'};
            ModelGRDusp100nM.stoichiometry = [-1,1,0,0,0,0,0,0,0,0,0,0,0,0;...
                1,-1,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,-1,1,1,-1,0,0,0,0,0,0,0,0;...
                0,0,1,-1,0,0,-1,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,1,-1,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,-1,1,1,-1,0;...
                0,0,0,0,0,0,0,0,0,1,-1,0,0,-1];
            ModelGRDusp100nM.useHybrid = true;
            ModelGRDusp100nM.hybridOptions.upstreamODEs = {'cytGR_a','nucGR_a','cytGR_b','nucGR_b'};
            ModelGRDusp100nM.solutionScheme = 'FSP';
            ModelGRDusp100nM.customConstraintFuns = [];
            ModelGRDusp100nM.fspOptions.bounds = [0;0;0;2;2;2;2;400];
            ModelGRDusp100nM.fittingOptions.modelVarsToFit = 1:4;
            ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('EricModDusp1');
            duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(1:4)).^2./(2*log10PriorStd(1:4).^2));
            ModelGRDusp100nM.fittingOptions.logPrior = duspLogPrior;
        
        %% STEP 2.A.2. -- Load pre-fit parameters into model.
        if loadPrevious
            load('EricModelDusp1_MMDex','DUSP1pars')
        else
            %% Pull the DUSP1 parameters from the Model
            % Find the indices of the desired parameter names
            knc_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'knc'));
            kg1_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'kg1'));
            gg1_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'gg1'));
            gg2_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'gg2'));
    
            % Extract the values and store them in DUSP1pars
            DUSP1pars = [ModelGRDusp100nM.parameters{knc_idx, 2}, ...
                        ModelGRDusp100nM.parameters{kg1_idx, 2}, ...
                        ModelGRDusp100nM.parameters{gg1_idx, 2}, ...
                        ModelGRDusp100nM.parameters{gg2_idx, 2}]; 
        end
        ModelGRDusp100nM.parameters(ModelGRDusp100nM.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        %% STEP 2.A.3. -- Load and Associate with DUSP1 smFISH Data (100nM Dex Only)
        % The commented code below would be needed to fit multiple conditions,
        % but that is not used in this case.  It is left here in case it is
        % needed in later stages of the project.
    
        if ~loadPrevious
            Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
            ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
            ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
            for i = 1:size(Dusp1FitCases,1)
                 ModelDusp1Fit{i} = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
                     {'rna','totalNucRNA'},...
                     {'Dex_Conc','100'});
                 ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
                 ModelDusp1parameterMap{i} = (1:4);
                 % Set Dex concentration.
                 ModelDusp1Fit{i}.parameters{13,2} = str2num(Dusp1FitCases{i,1});
                 ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
            end
            DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
        end
    
        ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
            {'rna','totalNucRNA'},{'Dex_Conc','100'});
        DUSP1pars = [ModelGRDusp100nM.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
    
        %%    STEP 2.B. -- Fit DUSP1 model at 100nM Dex.  
        % This will likely need to be rerun after Metropolis-Hastings Search (STEP 2.D.3).
        
        % Use fitMHiters to run until satisfied. TODO: Automate by computing statistics.
        fitMHiters = 2;
        
        for DS = 1:fitMHiters
            for i = 1:fitIters
                fitOptions.suppressFSPExpansion = true;
                DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
                    DUSP1pars, fitOptions);
                ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
                save('EricModelDusp1_MMDex','GRpars','DUSP1pars') 
            end
        
            %%    STEP 2.C. -- Plot predictions for other Dex concentrations.
            showCases = [1,1,1,1];
            makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1pars,Dusp1FitCases,showCases)
        
            %% STEP 2.D. -- Sample uncertainty for Dusp1 Parameters
            %   STEP 2.D.1. -- Compute sensitivity of the FSP solution
            ModelGRDusp100nM.solutionScheme = 'fspSens';
            sensSoln = ModelGRDusp100nM.solve();
            ModelGRDusp100nM.solutionScheme = 'FSP';
            %   STEP 2.D.2. -- Compute FIM
            %       define which species in model are not observed.
            ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};
            % compute the FIM
            fimResults = ModelGRDusp100nM.computeFIM(sensSoln.sens,'log');
            % In the following, the log-prior is used as a prior co-variance matrix.
            % This will be used in the FIM calculation as an FIM without new evidence 
            % being set equal to the inverse of this covariance matrix.  More rigorous
            % justification is needed to support this heuristic.
            fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
                diag(log10PriorStd(1:13).^2));
            FIMfree = fimTotal{1}(1:4,1:4);
            if min(eig(FIMfree))<1
                disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
                FIMfree = FIMfree + 1*eye(length(FIMfree));
            end
            covFree = FIMfree^-1;
            covFree = 0.5*(covFree+covFree');
        
            %%      STEP 2.D.3. -- Run Metropolis Hastings Search
            if loadPrevious
                MHDusp1File = 'MHDusp1_Dec92024';
                load(MHDusp1File)
            else
                MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
                MHFitOptions.thin=1;
                MHFitOptions.numberOfSamples=2000;
                MHFitOptions.burnIn=0;
                MHFitOptions.progress=true;
                MHFitOptions.numChains = 1;
                MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
                [DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
                    [], MHFitOptions, 'MetropolisHastings');
                delete('TMPEricMHDusp1.mat')
                ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
            end
            
            save('workspace_Feb4_2024.mat','ModelGRDusp100nM','DUSP1pars','fimTotal','sensSoln','combinedGRModel','MHResultsDusp1')
        
            %%      STEP 2.D.4. -- Plot the MH results
            figNew = figure;
            ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
            for j = 1:3
                for k = i:3
                    subplot(3,3,(j-1)*3+k)
                    CH = get(gca,'Children');
                    CH(1).Color=[1,0,1]; %
                    CH(1).LineWidth = 3;
                end
            end
        end
        
        %% Save results
        varNames = unique({'ModelGR'
            'GRfitCases'
            'log10PriorMean'
            'log10PriorStd'
            'GRpars'
            'MHResultsGR'
            'ModelGRparameterMap'
            'ModelGRfit'
            'boundGuesses'
            'ModelGRDusp100nM'
            'GRfitCases'
            'log10PriorMean'
            'log10PriorStd'
            'duspLogPrior'
            'DUSP1pars'
            'ModelGRfit'
            'fimResults'
            'MHResultsDusp1'
            'sensSoln'
            });
        
        save('workspaceDec9_2024',varNames{:}) % WARNING: THIS OVERWRITE THE PREVIOUSLY SAVED WORKSPACE - TODO: FIX
        case 1
        %% STEP 2.A.1. -- Add DUSP1 to the existing GR model.
            % Copy parameters from the 100nM Dex stim case in GR.
            fitIters = 3;
            ModelGRDusp100nM = ModelGRfit_a{3};
            ModelGRDusp100nM.species = {'offGene';'onGene';'cytGR_a';'nucGR_a';
                'rna';'cytGR_b';'nucGR_b'};
            ModelGRDusp100nM.initialCondition = [2;0;24;1;5;12;12];
            ModelGRDusp100nM.propensityFunctions = {'kon*offGene*nucGR_a';'koff*onGene*nucGR_b';
                '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex))*cytGR_a';'knc*nucGR_a';
                'kg1';'gg1*cytGR_a';'gg2*nucGR_a';'kr*onGene';'gr*rna';
                'kcn2*cyt_b';'knc*nucGR_b';'kg2';'gg1*cytGR_b';'gg2*nucGR_b'};
            ModelGRDusp100nM.stoichiometry = [-1,1,0,0,0,0,0,0,0,0,0,0,0,0;...
                1,-1,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,-1,1,1,-1,0,0,0,0,0,0,0,0;...
                0,0,1,-1,0,0,-1,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,1,-1,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,-1,1,1,-1,0;...
                0,0,0,0,0,0,0,0,0,1,-1,0,0,-1];
            ModelGRDusp100nM.useHybrid = true;
            ModelGRDusp100nM.hybridOptions.upstreamODEs = {'cytGR_a','nucGR_a','cytGR_b','nucGR_b'};
            ModelGRDusp100nM.solutionScheme = 'FSP';
            ModelGRDusp100nM.customConstraintFuns = [];
            ModelGRDusp100nM.fspOptions.bounds = [0;0;0;2;2;2;2;400];
            ModelGRDusp100nM.fittingOptions.modelVarsToFit = 1:4;
            ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('EricModDusp1');
            duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(1:4)).^2./(2*log10PriorStd(1:4).^2));
            ModelGRDusp100nM.fittingOptions.logPrior = duspLogPrior;
        
        %% STEP 2.A.2. -- Load pre-fit parameters into model.
        if loadPrevious
            load('EricModelDusp1_MMDex','DUSP1pars')
        else
            %% Pull the DUSP1 parameters from the Model
            % Find the indices of the desired parameter names
            knc_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'knc'));
            kg1_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'kg1'));
            gg1_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'gg1'));
            gg2_idx = find(strcmp(ModelGRDusp100nM.parameters(:,1), 'gg2'));
    
            % Extract the values and store them in DUSP1pars
            DUSP1pars = [ModelGRDusp100nM.parameters{knc_idx, 2}, ...
                        ModelGRDusp100nM.parameters{kg1_idx, 2}, ...
                        ModelGRDusp100nM.parameters{gg1_idx, 2}, ...
                        ModelGRDusp100nM.parameters{gg2_idx, 2}]; 
        end
        ModelGRDusp100nM.parameters(ModelGRDusp100nM.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        %% STEP 2.A.3. -- Load and Associate with DUSP1 smFISH Data (100nM Dex Only)
        % The commented code below would be needed to fit multiple conditions,
        % but that is not used in this case.  It is left here in case it is
        % needed in later stages of the project.
    
        if ~loadPrevious
            Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
            ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
            ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
            for i = 1:size(Dusp1FitCases,1)
                 ModelDusp1Fit{i} = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
                     {'rna','totalNucRNA'},...
                     {'Dex_Conc','100'});
                 ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
                 ModelDusp1parameterMap{i} = (1:4);
                 % Set Dex concentration.
                 ModelDusp1Fit{i}.parameters{13,2} = str2num(Dusp1FitCases{i,1});
                 ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
            end
            DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
        end
    
        ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
            {'rna','totalNucRNA'},{'Dex_Conc','100'});
        DUSP1pars = [ModelGRDusp100nM.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
    
        %%    STEP 2.B. -- Fit DUSP1 model at 100nM Dex.  
        % This will likely need to be rerun after Metropolis-Hastings Search (STEP 2.D.3).
        
        % Use fitMHiters to run until satisfied. TODO: Automate by computing statistics.
        fitMHiters = 2;
        
        for DS = 1:fitMHiters
            for i = 1:fitIters
                fitOptions.suppressFSPExpansion = true;
                DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
                    DUSP1pars, fitOptions);
                ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
                save('EricModelDusp1_MMDex','GRpars','DUSP1pars') 
            end
        
            %%    STEP 2.C. -- Plot predictions for other Dex concentrations.
            showCases = [1,1,1,1];
            makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1pars,Dusp1FitCases,showCases)
        
            %% STEP 2.D. -- Sample uncertainty for Dusp1 Parameters
            %   STEP 2.D.1. -- Compute sensitivity of the FSP solution
            ModelGRDusp100nM.solutionScheme = 'fspSens';
            sensSoln = ModelGRDusp100nM.solve();
            ModelGRDusp100nM.solutionScheme = 'FSP';
            %   STEP 2.D.2. -- Compute FIM
            %       define which species in model are not observed.
            ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};
            % compute the FIM
            fimResults = ModelGRDusp100nM.computeFIM(sensSoln.sens,'log');
            % In the following, the log-prior is used as a prior co-variance matrix.
            % This will be used in the FIM calculation as an FIM without new evidence 
            % being set equal to the inverse of this covariance matrix.  More rigorous
            % justification is needed to support this heuristic.
            fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
                diag(log10PriorStd(1:13).^2));
            FIMfree = fimTotal{1}(1:4,1:4);
            if min(eig(FIMfree))<1
                disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
                FIMfree = FIMfree + 1*eye(length(FIMfree));
            end
            covFree = FIMfree^-1;
            covFree = 0.5*(covFree+covFree');
        
            %%      STEP 2.D.3. -- Run Metropolis Hastings Search
            if loadPrevious
                MHDusp1File = 'MHDusp1_Dec92024';
                load(MHDusp1File)
            else
                MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
                MHFitOptions.thin=1;
                MHFitOptions.numberOfSamples=2000;
                MHFitOptions.burnIn=0;
                MHFitOptions.progress=true;
                MHFitOptions.numChains = 1;
                MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
                [DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
                    [], MHFitOptions, 'MetropolisHastings');
                delete('TMPEricMHDusp1.mat')
                ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
            end
            
            save('workspace_Feb4_2024.mat','ModelGRDusp100nM','DUSP1pars','fimTotal','sensSoln','combinedGRModel','MHResultsDusp1')
        
            %%      STEP 2.D.4. -- Plot the MH results
            figNew = figure;
            ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
            for j = 1:3
                for k = i:3
                    subplot(3,3,(j-1)*3+k)
                    CH = get(gca,'Children');
                    CH(1).Color=[1,0,1]; %
                    CH(1).LineWidth = 3;
                end
            end
        end
        
        %% Save results
        varNames = unique({'ModelGR'
            'GRfitCases'
            'log10PriorMean'
            'log10PriorStd'
            'GRpars'
            'MHResultsGR'
            'ModelGRparameterMap'
            'ModelGRfit'
            'boundGuesses'
            'ModelGRDusp100nM'
            'GRfitCases'
            'log10PriorMean'
            'log10PriorStd'
            'duspLogPrior'
            'DUSP1pars'
            'ModelGRfit'
            'fimResults'
            'MHResultsDusp1'
            'sensSoln'
            });
        
        save('workspaceDec9_2024',varNames{:}) % WARNING: THIS OVERWRITE THE PREVIOUSLY SAVED WORKSPACE - TODO: FIX
    end




    save('conv2solnTensor_postData','GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln', 'fspSolnsSMM',...
        'fimGR_a_withPrior','fimGR_b_withPrior','ModelGR_b_fimResults','fimGR_a_covFree','fimGR_b_covFree','MHResultsGR_a','MHFitOptions','ModelGRfit');
    
    

    %% compute likelihood
    %    logLikelihood = sum(log(conv2solnTensor{f}(PAs>0)).*PAs((PAs>0)),"all");
    %logLikelihood = sum(log(conv2solnTensor_postData{t}(fspsoln_sptensor_a_postData{t}>0)).*fspsoln_sptensor_a_postData{t}((fspsoln_sptensor_a_postData{t}>0)),"all");

    case 2
        %% Solve ODEs for fancy GR models
        disp('You have chosen GR-alpha + GR-beta by ODEs.')
        
        % Load experimental data
        %data = readtable("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv");  
        data = readtable("~/Selected_Columns_Reordered_Updated.csv"); % test 
        
        % Extract unique cell IDs and Dex concentrations
        cell_ids = unique(data.Cell_id);
        
        % Initialize storage for simulation results
        all_simulations = table();
        
        %for i = 1:length(cell_ids)
        for i = 1:10 % WARNING:  hacky sack!
            % Extract data for the current cell
            cell_data = data(data.Cell_id == cell_ids(i), :);
        
            % Extract time and measurements
            time_data = cell_data.Time_index;
            normgrnuc = cell_data.normgrnuc; % nucGR_a + nucGR_b
            normgrcyt = cell_data.normgrcyt; % cytGR_a + cytGR_b
        
            % Initial conditions for this cell
            cytGR_a_0 = 5.0; % Initial concentration of cytGR_a
            nucGR_a_0 = 1.0; % Initial concentration of nucGR_a
            cytGR_b_0 = 5.0; % Initial concentration of cytGR_b
            nucGR_b_0 = 1.0; % Initial concentration of nucGR_b
            y0 = [cytGR_a_0; nucGR_a_0; cytGR_b_0; nucGR_b_0];
        
            % Time span for this cell
            tspan = [min(time_data), max(time_data)];
        
            % Solve the ODEs
            [t, y] = ode45(@crs_odes, tspan, y0);
        
            % Extract results
            cytGR_a = y(:, 1); % cytGR_a
            nucGR_a = y(:, 2); % nucGR_a
            cytGR_b = y(:, 3); % cytGR_b
            nucGR_b = y(:, 4); % nucGR_b
        
            % Compute aggregate values for comparison
            sim_normgrnuc = nucGR_a + nucGR_b; % nucGR_a + nucGR_b
            sim_normgrcyt = cytGR_a + cytGR_b; % cytGR_a + cytGR_b
        
            % Store results in a table
            sim_table = table(t, sim_normgrnuc, sim_normgrcyt, repmat(cell_ids(i), length(t), 1), ...
                              'VariableNames', {'time', 'sim_normgrnuc', 'sim_normgrcyt', 'Cell_id'});
            
            % Append to all simulations
            all_simulations = [all_simulations; sim_table];
        
            % Plot results for the current cell
            figure;
            plot(t, sim_normgrnuc, 'LineWidth', 2, 'DisplayName', 'Simulated normgrnuc');
            hold on;
            plot(t, sim_normgrcyt, 'LineWidth', 2, 'DisplayName', 'Simulated normgrcyt');
            plot(time_data, normgrnuc, 'o', 'DisplayName', 'Experimental normgrnuc');
            plot(time_data, normgrcyt, 'o', 'DisplayName', 'Experimental normgrcyt');
            hold off;
        
            title(['Comparison for Cell ID: ', num2str(cell_ids(i))]);
            xlabel('Time');
            ylabel('Concentration');
            legend show;
            grid on;
        end
end

%% ODE system for GR_alpha + GR_beta
% Define the ODE system
function dydt = crs_odes(t, y)

    % Define parameters
    kcn0 = 0.005; % base rate of cytGR_a -> nucGR_a
    kcn1 = 0.02; % rate of cytGR_a -> nucGR_a with Dex
    MDex = 5.0;
    Dex0 = 100.0; %Get Dex concentration from input file
    gDex = 0.003;
    
    ba1 = 14e-5; % birth of cytGR_a
    da1 = 1e-5; % death of cytGR_a
    ka1 = 0.01; % nucGR_a -> cytGR_a
    da2 = 1e-6; % death of nucGR_a
    
    bb1 = 14e-5; % birth of cytGR_b
    db1 = 1e-5; % death of cytGR_b
    kb1 = 0.01; % cytGR_b -> nucGR_b
    kb2 = 0.01; % nucGR_b -> cytGR_b
    db2 = 1e-6; % death of nucGR_b
    
    % Define input signal
    IDex = @(t) Dex0 * exp(-gDex * t);

    % Unpack variables
    cytGR_a = y(1); % cytGR_a
    nucGR_a = y(2); % nucGR_a
    cytGR_b = y(3); % cytGR_b
    nucGR_b = y(4); % nucGR_b

    % Reaction rates
    w1 = (kcn0 + kcn1 * IDex(t) / (MDex + IDex(t))) * cytGR_a;
    w2 = ba1;
    w3 = da1 * cytGR_a;
    w4 = ka1 * nucGR_a;
    w5 = da2 * nucGR_a;

    w6 = bb1;
    w7 = db1 * cytGR_b;
    w8 = kb1 * cytGR_b;
    w9 = kb2 * nucGR_b;
    w10 = db2 * nucGR_b;

    % ODEs
    dcytGRa_dt = -w1 + w2 - w3 + w4;
    dnucGR_a_dt = w1 - w4 - w5;
    dcytGR_b_dt = w6 - w7 - w8 + w9;
    dnucGR_b_dt = w8 - w9 - w10;

    % Return derivatives
    dydt = [dcytGRa_dt; dnucGR_a_dt; dcytGR_b_dt; dnucGR_b_dt];
end



%% ODE system for GR_alpha + GR_beta
% Define the ODE system
function dydt = gragrb(t, y)

    % Define parameters
    kcn0 = 0.005; % base rate of cytGR_a -> nucGR_a
    kcn1 = 0.02; % rate of cytGR_a -> nucGR_a with Dex
    MDex = 5.0;
    Dex0 = 100.0; %Get Dex concentration from input file
    gDex = 0.003;
    
    ba1 = 14e-5; % birth of cytGR_a
    da1 = 1e-5; % death of cytGR_a
    ka1 = 0.01; % nucGR_a -> cytGR_a
    da2 = 1e-6; % death of nucGR_a
    
    bb1 = 14e-5; % birth of cytGR_b
    db1 = 1e-5; % death of cytGR_b
    kb1 = 0.01; % cytGR_b -> nucGR_b
    kb2 = 0.01; % nucGR_b -> cytGR_b
    db2 = 1e-6; % death of nucGR_b
    
    % Define input signal
    IDex = @(t) Dex0 * exp(-gDex * t);

    % Unpack variables
    cytGR_a = y(1); % cytGR_a
    nucGR_a = y(2); % nucGR_a
    cytGR_b = y(3); % cytGR_b
    nucGR_b = y(4); % nucGR_b

    % Reaction rates
    w1 = (kcn0 + kcn1 * IDex(t) / (MDex + IDex(t))) * cytGR_a;
    w2 = ba1;
    w3 = da1 * cytGR_a;
    w4 = ka1 * nucGR_a;
    w5 = da2 * nucGR_a;

    w6 = bb1;
    w7 = db1 * cytGR_b;
    w8 = kb1 * cytGR_b;
    w9 = kb2 * nucGR_b;
    w10 = db2 * nucGR_b;

    % ODEs
    dcytGRa_dt = -w1 + w2 - w3 + w4;
    dnucGR_a_dt = w1 - w4 - w5;
    dcytGR_b_dt = w6 - w7 - w8 + w9;
    dnucGR_b_dt = w8 - w9 - w10;

    % Return derivatives
    dydt = [dcytGRa_dt; dnucGR_a_dt; dcytGR_b_dt; dnucGR_b_dt];
end

%% ODE system for GR_alpha + GR_beta
% Define the ODE system
function dydt = gr_a_gr_b(t, y)

    % Define parameters
    kcn0 = 0.005; % base rate of cytGR_a -> nucGR_a
    kcn1 = 0.02; % rate of cytGR_a -> nucGR_a with Dex
    MDex = 5.0;
    Dex0 = 100.0; %Get Dex concentration from input file
    gDex = 0.003;
    
    ba1 = 14e-5; % birth of cytGR_a
    da1 = 1e-5; % death of cytGR_a
    ka1 = 0.01; % nucGR_a -> cytGR_a
    da2 = 1e-6; % death of nucGR_a
    
    bb1 = 14e-5; % birth of cytGR_b
    db1 = 1e-5; % death of cytGR_b
    kb1 = 0.01; % cytGR_b -> nucGR_b
    kb2 = 0.01; % nucGR_b -> cytGR_b
    db2 = 1e-6; % death of nucGR_b
    
    % Define input signal
    IDex = @(t) Dex0 * exp(-gDex * t);

    % Unpack variables
    cytGR_a = y(1); % cytGR_a
    nucGR_a = y(2); % nucGR_a
    cytGR_b = y(3); % cytGR_b
    nucGR_b = y(4); % nucGR_b

    % Reaction rates
    w1 = (kcn0 + kcn1 * IDex(t) / (MDex + IDex(t))) * cytGR_a;
    w2 = ba1;
    w3 = da1 * cytGR_a;
    w4 = ka1 * nucGR_a;
    w5 = da2 * nucGR_a;

    w6 = bb1;
    w7 = db1 * cytGR_b;
    w8 = kb1 * cytGR_b;
    w9 = kb2 * nucGR_b;
    w10 = db2 * nucGR_b;

    % ODEs
    dcytGRa_dt = -w1 + w2 - w3 + w4;
    dnucGR_a_dt = w1 - w4 - w5;
    dcytGR_b_dt = w6 - w7 - w8 + w9;
    dnucGR_b_dt = w8 - w9 - w10;

    % Return derivatives
    dydt = [dcytGRa_dt; dnucGR_a_dt; dcytGR_b_dt; dnucGR_b_dt];
end
