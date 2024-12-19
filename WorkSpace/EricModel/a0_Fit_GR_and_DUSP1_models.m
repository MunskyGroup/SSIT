%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations
close all 
clear
addpath(genpath('../../src'));
loadPrevious = false;
savedWorkspace = 'workspaceDec9_2024';
addpath('tmpPropensityFunctions');

%% GR model
%% STEP 0 -- Preliminaries: Load previously computed data or define GR model from scratch
if loadPrevious
    % STEP 0.A. -- Option to load previous workspace from file -- just for testing.
    varNames = {'ModelGR'
    'GRfitCases'
    'log10PriorMean'
    'log10PriorStd'
    'GRpars'
    'ModelGRparameterMap'
    'ModelGRfit'
    'boundGuesses'
    'ModelGRDusp100nM'
    'GRfitCases'
    'log10PriorMean'
    'log10PriorStd'
    'duspLogPrior'
    'DUSP1pars'
    'Dusp1FitCases'
    'ModelGRfit'
    'extendedMod'
    'ModelGRDusp100nM_ext_red'
    'fullPars'
    'fimResults'
    'MHResultsGR'
    'MHResultsDusp1'
    'sensSoln.sens'
    };
    load(savedWorkspace,varNames{:})
    fitOptions = optimset('Display','iter','MaxIter',10);
    try
        ModelGRfit{1}.propensitiesGeneral{1}.stateDependentFactor(0);
    catch
        for i=1:length(ModelGRfit)
            ModelGRfit{i} = ModelGRfit{i}.formPropensitiesGeneral('GR_Model');
        end
    end
    fitIters = 10;
    load('EricModel_MMDex','GRpars')
    ModelGR.parameters(5:12,2) = num2cell([GRpars]);
    GRpars
    %GRpars = GRpars';
else
    %% STEP 0.B. -- Create Base Model for GR Only
    %    Here, I set up the model for the GR translocation dynamics.
    fitOptions = optimset('Display','iter','MaxIter',300);

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
end

%% STEP 1 -- Fit GR Models.  
% STEP 1 will need to be rerun until satisfied.  Use fitMHiters as needed.
% TODO: Automate with statistics.
% Set for STEP1 -- Fit GR Models
fitIters = 30;
fitMHiters = 20;

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
save('workspaceDec9_2024.mat','GRpars', 'ModelGRfit', 'combinedGRModel','MHResultsGR', 'log10PriorStd')

%%  STEP 2 -- Extend Model to Include DUSP1 Activation, Production, and Degradation
if loadPrevious
    varNamesDUSP1 = {'ModelGR'
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
    'MHResultsDusp1'};
    load(savedWorkspace,varNamesDUSP1{:})
    fitOptions = optimset('Display','iter','MaxIter',10);
    fitIters = 1;
    try
        ModelGRDusp100nM.propensitiesGeneral{1}.stateDependentFactor(0);
    catch
        ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('DUSP1_Model');
    end
else
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
end
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

save('workspaceDec9_2024',varNames{:}) % WARNING: THIS OVERWRITE THE PREVIOUSLY SAVED WORKSPACE - TODO: FIX

%%  STEP 3. -- Model Extensions using ODE Analyses
if loadPrevious
    vaNamesExtended = {'ModelGRfit'
        'extendedMod'
        'ModelGRDusp100nM_ext_red'
        'fullPars'
        };
    load(savedWorkspace,vaNamesExtended{:})
    try
        extendedMod.propensitiesGeneral{1}.stateDependentFactor(0);
    catch
        extendedMod = extendedMod.formPropensitiesGeneral('NucCyt_Model');
    end
    try
        ModelGRDusp100nM_ext_red.propensitiesGeneral{1}.stateDependentFactor(0);
    catch
        ModelGRDusp100nM_ext_red = ModelGRDusp100nM_ext_red.formPropensitiesGeneral('ExtModel100nm');
    end
    try
        for i=1:length(ModelGRfit)
            ModelGRfit{i}.propensitiesGeneral{1}.stateDependentFactor(0);
        end
    catch
        for i=1:length(ModelGRfit)
            ModelGRfit{i} = ModelGRfit{i}.formPropensitiesGeneral(['GRModB_',num2str(i)]);
        end
    end
else
    %%    STEP 3.A.1 -- Extend model to include nuclear and cytoplasmic RNA
    extendedMod = ModelGRDusp100nM;
    extendedMod = extendedMod.addSpecies({'rCyt'},0);
    % Adjust the final reaction to be nuc -> cyt transport.
    extendedMod.propensityFunctions{9,1} = 'knuc2cyt*rna';
    extendedMod.parameters{4,1} = 'knuc2cyt';
    extendedMod.stoichiometry(:,9) = 0;
    extendedMod.stoichiometry([5,6],9) = [-1;1];
    % Add nuc degradation reaction.
    extendedMod.propensityFunctions{10,1} = 'knucdeg*rna';
    extendedMod.stoichiometry(:,10) = 0;
    extendedMod.stoichiometry(5,10) = -1;
    extendedMod.parameters(14,:) = {'knucdeg',0.001};
    % Add cytoplasmic degradation as a reaction
    extendedMod.propensityFunctions{11,1} = 'degCyt*rCyt';
    extendedMod.stoichiometry(:,11) = 0;
    extendedMod.stoichiometry(6,11) = -1;
    extendedMod.parameters(15,:) = {'degCyt',100};
    % The following code can be used to verify that the model changes so far
    % are consistent with the previous FSP model.
    % This is useful to verify before changing the reactions too much.
    % However, we do not expect the FSP analysis to be very fast just yet, so
    % it would not make much sense to run this for fitting.
    if false
        extendedMod.initialCondition(6) = 1;
        extendedMod.fspOptions.bounds = [0;0;0;0;2;2;400;1]
        extendedMod.fspOptions.verbose = 1;
        extendedMod = extendedMod.formPropensitiesGeneral('extendedModel');
        [~,extendedMod.fspOptions.bounds] = extendedMod.solve;
        extendedMod = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
            {'rna','RNA_DUSP1_nuc'; ...
            'rCyt','RNA_DUSP1_cyto'},...
            {'Dex_Conc','100'});
        extendedMod.makeFitPlot
    end
    %%    STEP 3.A.2 -- Create reduced model for Nuclear DUSP1 but with same parameters as new model
    ModelGRDusp100nM_ext_red = ModelGRDusp100nM;
    % Adjust the final reaction to be nuc -> cyt transport.
    ModelGRDusp100nM_ext_red.propensityFunctions{9,1} = 'knuc2cyt*rna';
    ModelGRDusp100nM_ext_red.parameters{4,1} = 'knuc2cyt';
    ModelGRDusp100nM_ext_red.stoichiometry(:,9) = 0;
    ModelGRDusp100nM_ext_red.stoichiometry(5,9) = -1;
    % Add nuc degradation reaction.
    ModelGRDusp100nM_ext_red.propensityFunctions{10,1} = 'knucdeg*rna';
    ModelGRDusp100nM_ext_red.stoichiometry(:,10) = 0;
    ModelGRDusp100nM_ext_red.stoichiometry(5,10) = -1;
    ModelGRDusp100nM_ext_red.parameters(14,:) = {'knucdeg',0.001};
    ModelGRDusp100nM_ext_red = ModelGRDusp100nM_ext_red.formPropensitiesGeneral('ModelGRDusp100nM_ext_red');
    ModelGRDusp100nM_ext_red.fittingOptions.modelVarsToFit = [1:4,14];

    %log10PriorMean = [-1 -1 0 -2,...
    %    -1 -3 -2 -1 -2 -2 -2 0.5, 2];
log10PriorMean = [-1 -1 0 -2,... %dusp1 pars
    -1 -3 -2 -1 -2 -2 -2 0.5, ...%GR pars
    2, ... % Dex concentration -- known
    -2, -3]; % dusp1 transport, cyt RNA degradation
log10PriorStd = 2*ones(1,14);
ModelGRDusp100nM_ext_red.fittingOptions.logPrior = @(x)-sum((log10(x)-log10PriorMean([1:4,14])).^2./(2*log10PriorStd([1:4,14]).^2));
    
    ModelGRDusp100nM_ext_red.makeFitPlot

    %%    STEP 3.B.1 -- Load data into the model
    extendedMod = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','RNA_DUSP1_nuc'; ...
        'rCyt','RNA_DUSP1_cyto'},...
        {'Dex_Conc','100'});
    %%    STEP 3.B.2 -- Switch solver to ODE and generate model codes
    extendedMod.solutionScheme = 'ode';
    extendedMod.useHybrid = false;
    extendedMod = extendedMod.formPropensitiesGeneral('ODEEricCyt');
    %%    STEP 3.C. -- Change parameters manually, solve, and make plots
    %extendedMod.parameters(4,:) = {'knuc2cyt',0.027};
    %extendedMod.parameters(14,:) = {'knucdeg',0.001};
    %extendedMod.parameters(15,:) = {'degCyt',0.01};
    extendedMod.parameters(4,:) = {'knuc2cyt',0.06};
    extendedMod.parameters(14,:) = {'knucdeg',0.005}; % TODO: check this... elsewhere in log10PriorMean it says the 14th parameter is dusp1 transport.....
    extendedMod.parameters(15,:) = {'degCyt',0.05};
    soln = extendedMod.solve;
    plotODEresults(extendedMod,soln,ModelGRfit{3})
    %%    STEP 3.D.1. -- Fit new parameters to match all ODE data.
    extendedMod.fittingOptions.modelVarsToFit = [4,14,15];
    extendedMod.fittingOptions.logPrior = [];
    pars = [extendedMod.parameters{extendedMod.fittingOptions.modelVarsToFit,2}];
    pars = extendedMod.maximizeLikelihood(pars,fitOptions);
    %%    STEP 3.D.2. -- Plot ODE fit results.
    extendedMod.parameters(extendedMod.fittingOptions.modelVarsToFit,2) = num2cell(pars);
    soln = extendedMod.solve;
    plotODEresults(extendedMod,soln,ModelGRfit{3})
    %%    STEP 3.F.1. -- Create New Objective Function Combining all of the previous ones.
    % Remove all priors from individual models.
    ModelGRDusp100nM_ext_red.fittingOptions.logPrior = [];
    for i=1:3
        ModelGRfit{1}.fittingOptions.logPrior = [];
        ModelGRfit{1}.tSpan = ModelGRfit{1}.dataSet.times;
    end
    extendedMod.fittingOptions.logPrior = [];
    fullPars = [extendedMod.parameters{[1:12,14,15],2}];
    % otherewise, we will use the set that was saved in the data dump from
    % the older workspace.
end
%%    STEP 3.F.2. -- Fit all objective functions at once.
% Create prior for all parameters
log10PriorMean = [-1 -1 0 -2,... %dusp1 pars
    -1 -3 -2 -1 -2 -2 -2 0.5, ...%GR pars
    NaN, ... % Dex concentration -- known
    -2, -3]; % dusp1 transport, cyt RNA degradation
log10PriorStd = 2*ones(1,15);
logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:12,14,15])).^2./(2*log10PriorStd([1:12,14,15]).^2));
% extendedMod.fittingOptions.modelVarsToFit = [1:12,14,15];
Organization = {ModelGRfit{1},[5:12],[5:12],'computeLikelihood',1;...
    ModelGRfit{2},[5:12],[5:12],'computeLikelihood',1;...
    ModelGRfit{3},[5:12],[5:12],'computeLikelihood',1;...
    ModelGRDusp100nM_ext_red,[1:12,14],[1:13],'computeLikelihood',1;...
    extendedMod,[1:12,14:15],[1:14],'computeLikelihoodODE',0.01};
Organization = getTotalFitErr(Organization,fullPars,true);
getTotalFitErr(Organization,fullPars,false)
objAll = @(x)-getTotalFitErr(Organization,exp(x),false)-logPriorAll(exp(x));
for jj = 1:fitIters
    fullParsLog = log(fullPars);
    fullPars = exp(fminsearch(objAll,fullParsLog,fitOptions));
end

%%    STEP 3.F.3. -- Plot all results
% Plot GR Distribution
makeGRPlots(combinedGRModel,fullPars(5:12))

% Plot DUSP1 100nm FIT and other PREDICTED distributions
showCases = [1,1,1,1];
ModelGRDusp100nM_ext_red.fittingOptions.modelVarsToFit = [1:12,14];
ModelGRDusp100nM_ext_red.parameters([1:12,14],2) = num2cell(fullPars([1:13]));
makePlotsDUSP1({ModelGRDusp100nM_ext_red},ModelGRDusp100nM_ext_red,fullPars([1:13]),Dusp1FitCases,showCases)

% Plot ODE Cyt FIT Results at 100nM Dex
extendedMod.parameters([1:12,14,15],2) = num2cell(fullPars);
soln100 = extendedMod.solve;
plotODEresults(extendedMod,soln100,ModelGRfit{3},500)
set(gcf,'Name','ODE Fits -- 100nM Dex')

% Plot ODE Predictions at other DEX concentrations
extendedMod0p3 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
    {'rna','RNA_DUSP1_nuc'; ...
    'rCyt','RNA_DUSP1_cyto'},...
    {'Dex_Conc','0.3'});
extendedMod0p3.parameters(13,:) = {'Dex0',0.3};

extendedMod1 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
    {'rna','RNA_DUSP1_nuc'; ...
    'rCyt','RNA_DUSP1_cyto'},...
    {'Dex_Conc','1'});
extendedMod1.parameters(13,:) = {'Dex0',1.0};

extendedMod10 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
    {'rna','RNA_DUSP1_nuc'; ...
    'rCyt','RNA_DUSP1_cyto'},...
    {'Dex_Conc','10'});
extendedMod10.parameters(13,:) = {'Dex0',10};

plotODEresults(extendedMod1,extendedMod1.solve,ModelGRfit{1},501)
set(gcf,'Name','ODE Predictions -- 1.0nM Dex')

plotODEresults(extendedMod10,extendedMod10.solve,ModelGRfit{2},502)
set(gcf,'Name','ODE Predictions -- 10nM Dex')

%%    STEP 3.G.1. -- Create SSA model to predict cytoplasmic distributions
SSAModel_100 = extendedMod;
SSAModel_100.solutionScheme = 'SSA';
SSAModel_100.tSpan = [-500,ModelGRDusp100nM.tSpan];
% A negative initial time is needed to allow model to equilibrate before
% starting.  This causes long run times.
SSAModel_100.initialCondition = [2;0;round(soln100.ode(1,3:6))'];
SSAModel_100.initialTime = SSAModel_100.tSpan(1);
SSAModel_100.useHybrid = false;
SSAModel_100.ssaOptions.useParalel = true;
SSAModel_100 = SSAModel_100.formPropensitiesGeneral('SSAExtendedModel');

%%    STEP 3.G.2. -- Run SSA Simulations
ssaSoln_100 = SSAModel_100.solve;

%%    STEP 3.G.3. -- Plot SSA Results for Cytoplasmic Distributions (100nM Dex)
% Plot Fits for 100nM Dex 
makeCytDistPlots(ssaSoln_100,extendedMod,600,[2:10],[1:9],6,2)

%%    STEP 3.F.1. -- Predict Cyt distributions for other 0.3nM Dex 
SSAModel_0p3 = SSAModel_100;
SSAModel_0p3.parameters(13,:) = {'Dex0',0.3};
SSAModel_0p3.tSpan = [-500,extendedMod0p3.tSpan];
ssaSoln_0p3 = SSAModel_0p3.solve;

%%    STEP 3.F.2. -- Predict Cyt distributions for other 1.0nM Dex 
SSAModel_1 = SSAModel_100;
SSAModel_1.parameters(13,:) = {'Dex0',1.0};
SSAModel_1.tSpan = [-500,extendedMod1.tSpan];
ssaSoln_1 = SSAModel_1.solve;

%%    STEP 3.F.3. -- Predict Cyt distributions for other 10nM Dex 
SSAModel_10 = SSAModel_100;
SSAModel_10.parameters(13,:) = {'Dex0',10};
SSAModel_10.tSpan = [-500,extendedMod10.tSpan];
ssaSoln_10 = SSAModel_10.solve;

%%    STEP 3.F.4. -- Make resulting plots
makeCytDistPlots(ssaSoln_0p3,extendedMod0p3,601,[2:8],[1:7],6,2)
makeCytDistPlots(ssaSoln_1,extendedMod1,602,[3:8],[1:6],6,2)
makeCytDistPlots(ssaSoln_10,extendedMod10,603,[3:8],[1:6],6,2)

%%
%% Save Results for Easier Use in subsequent runs.
%parsAll_GR_Dusp1_TS = [extendedMod.parameters{:,2}];
%parsAll_GR_Dusp1_TS(16) = parsAllandTS(end);
varNames = unique({'ModelGR'
    'GRfitCases'
    'log10PriorMean'
    'log10PriorStd'
    'GRpars'
    'ModelGRparameterMap'
    'ModelGRfit'
    'boundGuesses'
    'ModelGRDusp100nM'
    'GRfitCases'
    'log10PriorMean'
    'log10PriorStd'
    'duspLogPrior'
    'DUSP1pars'
    'Dusp1FitCases'
    'ModelGRfit'
    'extendedMod'
    'ModelGRDusp100nM_ext_red'
    'fullPars'
    'fimResults'
    'MHResultsGR'
    'MHResultsDusp1'
    'sensSoln'
    });

save('workspaceDec9_2024',varNames{:})

%% Extra Function

function makeCytDistPlots(ssaSoln_100,extendedMod,fignum,timeIndsMod,timeIndsDat,speciesIndMod,speciesIndDat)
arguments
    ssaSoln_100
    extendedMod
    fignum = 1;
    timeIndsMod = [];
    timeIndsDat = [];
    speciesIndMod = 1;
    speciesIndDat = 1;
end
figure(fignum); clf;
if isempty(timeIndsMod)
    timeIndsMod = [1:length(ssaSoln_100.T_array)];
end
if isempty(timeIndsDat)
    timeIndsDat = [1:length(extendedMod.dataSet.times)];
end

if length(timeIndsMod)~=length(timeIndsDat)
    error('Length of data and model time points must be equal')
end

Nrows = ceil(sqrt(length(timeIndsMod)));
for i = 1:length(timeIndsMod)
    subplot(Nrows,Nrows,i)

    % Add SSA to histogram plot
    M = squeeze(ssaSoln_100.trajs(speciesIndMod,timeIndsMod(i),:));
    H = histogram(M,'Normalization','pdf');
    hold on

    % Add data to histogram plot
    dMat = double(extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(timeIndsDat(i),:,:));
    N = sum(dMat,"all");
    if speciesIndDat==1
        PD = [0;sum(dMat,2)/N];
    else
        PD = [0;sum(dMat,1)'/N];
    end
    
    binEdges = round(H.BinEdges);
    nBins = length(binEdges)-1;
    PDbinned = zeros(nBins,1);
    binwidth = binEdges(2)-binEdges(1);
    for j = 1:nBins
        PDbinned(j) = sum(PD(binEdges(j)+1:binEdges(j+1)));
    end

    PDbinned = [PDbinned;1-sum(PDbinned)];
    stairs(binEdges,PDbinned/binwidth,'linewidth',2)

    set(gca,'FontSize',15,'ylim',[0,0.03])
    
end
end
