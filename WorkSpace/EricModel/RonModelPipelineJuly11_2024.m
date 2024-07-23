%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are exp cted that would result in
%   slightly different parameter combinations
close all 
clear all
addpath(genpath('../../src'));

loadPrevious = true;
savedWorkspace = 'workspaceJuly22';
addpath('tmpPropensityFunctions');

%% STEP 0 -- Preliminaries.
if loadPrevious
    %% STEP 0.A -- Option to load previous workspace from file -- just for testing.
    varNames = {'ModelGR',
    'GRfitCases'
    'log10PriorMean'
    'log10PriorStd'
    'GRpars'
    'ModelGRparameterMap'
    'ModelGRfit'
    'boundGuesses'};

    load(savedWorkspace,varNames{:})
    fitOptions = optimset('Display','iter','MaxIter',10);
    fitIters = 1; 

    try
        ModelGRfit{1}.propensitiesGeneral{1}.stateDependentFactor(0);
    catch
        for i=1:length(ModelGRfit)
            ModelGRfit{i} = ModelGRfit{i}.formPropensitiesGeneral('GR_Model');
        end
    end
else
    fitOptions = optimset('Display','iter','MaxIter',300);
    %% STEP 0.B.1. -- Create Base Model for GR Only
    % Here, I set up the model for the GR translocation dynamics.
    ModelGR = SSIT;
    ModelGR.species = {'cytGR';'nucGR'};
    ModelGR.initialCondition = [20;1];

    ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';...
        'kg1';'gg1*cytGR';'gg2*nucGR'};
    ModelGR.stoichiometry = [-1,1,1,-1,0;...
        1,-1,0,0,-1];

    ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
        'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',14e-5;...
        'gg1',1e-5;'gg2',1e-6;'MDex',5;'Dex0',100});

    log10PriorMean = [-1 -1 0 -2,...
        -1 -3 -2 -1 -2 -2 -2 0.5, 2];
    log10PriorStd = 2*ones(1,13);
    % the log prior will be applied to the fit to multiple models as an
    % additional constraint.
    
    ModelGR.fittingOptions.logPrior = [];  
    % So it is left out of the prior, since we only want it to be calculated once.

    ModelGR.fspOptions.initApproxSS = true;

    ModelGR.fittingOptions.modelVarsToFit = (5:12);
    
    ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
    ModelGR = ModelGR.formPropensitiesGeneral('EricRonModGR');
    ModelGR.customConstraintFuns = {'cytGR+nucGR'};
    % TODO - Alex - Constraints for removing stiff dimensions.
    
    [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve;
    [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve(FSPGrSoln.stateSpace);
    %% STEP 0.B.2. -- Load previously fit parameter values (optional)
    load('EricModel_MMDex','GRpars')
    ModelGR.parameters(5:12,2) = num2cell([GRpars]);

    %% STEP 0.B.3. -- Associate GR Data with Different Instances of Model (10,100nm Dex)
    GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
        '10','10',102,'GR Fit (10nM Dex)';...
        '100','100',103,'GR Fit (100nM Dex)'};

    ModelGRparameterMap = cell(1,size(GRfitCases,1));
    ModelGRfit = cell(1,size(GRfitCases,1));
    % ModelGRODEfit = cell(1,size(GRfitCases,1));
    for i=1:3
        ModelGRfit{i} = ModelGR.loadData("EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
            {'nucGR','normgrnuc';'cytGR','normgrcyt'},...
            {'Condition','GR_timesweep';'Dex_Conc',GRfitCases{i,2}});
        ModelGRfit{i}.parameters(13,:) = {'Dex0', str2num(GRfitCases{i,1})};
        ModelGRparameterMap(i) = {(1:8)};
        % parameters 1 - 8 refer to the parameter set that is relevant to
        % the entire class of models.  In this case, these refer to
        % global parameters 5:15 (GR parameters).
    end
    %% STEP 0.B.4. -- Make Guesses for the FSP bounds
    % This is sometimes necessary when using an uninduced steady state as the
    % initial condition. You need to guess a reasonalbe statespace or the
    % computation of the SS can be inaccurate.
    ModelGR.customConstraintFuns = {'cytGR+nucGR'};
    for i = 1:3
        boundGuesses{i} = [0;0;30;30;30];
        % First N are lower bounds.  Next N is upper bound.  Remaining are
        % custom.
    end
end
%%
%% STEP 1 -- GR Model
%%     STEP 1.A. -- Combine all three GR models and fit using a single parameter set.
for i = 1:3
    ModelGRfit{i}.tSpan = ModelGRfit{i}.dataSet.times;
end

logPriorGR = @(x)-sum((log10(x)-log10PriorMean(5:12)).^2./(2*log10PriorStd(5:12).^2));
for jj = 1:fitIters
    combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap,logPriorGR);
    combinedGRModel = combinedGRModel.initializeStateSpaces(boundGuesses);
    combinedGRModel = combinedGRModel.updateModels(GRpars,false);
    GRpars = combinedGRModel.maximizeLikelihood(...
        GRpars, fitOptions);
    save('EricModel_MMDex','GRpars') 
end
%%     STEP 1.B. -- Compute FIM for GR parameters.
combinedGRModel = combinedGRModel.computeFIMs([],'log');
fimGR_withPrior = combinedGRModel.FIM.totalFIM+... % the FIM in log space.
    diag(1./(log10PriorStd(ModelGR.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.

%%     STEP 1.C. -- Run MH on GR Models.
% MHFitOptions.thin=1;
% MHFitOptions.numberOfSamples=100;
% MHFitOptions.burnIn=0;
% MHFitOptions.progress=true;
% MHFitOptions.numChains = 1;
% MHFitOptions.useFIMforMetHast = true;
% MHFitOptions.saveFile = 'TMPEricMHGR.mat';
% [~,~,MHResultsGR] = combinedGRModel.maximizeLikelihood(...
%     GRpars, MHFitOptions, 'MetropolisHastings');
% delete(MHFitOptions.saveFile)
 
%%     STEP 1.D. -- Make Plots of GR Fit Results
makeGRPlots(combinedGRModel,GRpars)

%%
%%  STEP 2 -- Extend Model to Include DUSP1 Activation, Production, and Degradation
if loadPrevious
    varNamesDUSP1 = {'ModelGRDusp100nM'
    'GRfitCases'
    'log10PriorMean'
    'log10PriorStd'
    'duspLogPrior'
    'DUSP1pars'
    'Dusp1FitCases'};
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
    ModelGRDusp100nM.fspOptions.bounds = [0;0;0;2;2;400];
    ModelGRDusp100nM.fittingOptions.modelVarsToFit = 1:4;
    ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('EricModDusp1');
    duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(1:4)).^2./(2*log10PriorStd(1:4).^2));
    ModelGRDusp100nM.fittingOptions.logPrior = duspLogPrior;

    %% STEP 2.A.2. -- Load pre-fit parameters into model.
    load('EricModelDusp1_MMDex','DUSP1pars')
    ModelGRDusp100nM.parameters(ModelGRDusp100nM.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);

    %% STEP 2.A.3. -- Load and Associate with DUSP1 smFISH Data (100nM Dex Only)
    % The commented code below would be needed to fit multiple conditions,
    % but that is not used in this case.  It is left here in case it is
    % needed in later stages of the project.
    % % % Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
    % % % ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
    % % % ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
    % % % for i = 1:size(Dusp1FitCases,1)
    % % %     ModelDusp1Fit{i} = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
    % % %         {'rna','totalNucRNA'},...
    % % %         {'Dex_Conc','100'});
    % % %     ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
    % % % 
    % % %     ModelDusp1parameterMap{i} = (1:4);
    % % %     % Set Dex concentration.
    % % %     ModelDusp1Fit{i}.parameters{13,2} = str2num(Dusp1FitCases{i,1});
    % % %     ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
    % % % end
    % % % DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];

    ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','totalNucRNA'},{'Dex_Conc','100'});
    DUSP1pars = [ModelGRDusp100nM.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
end

%%    STEP 2.B. -- Fit DUSP1 model at 100nM Dex.
for i = 1:fitIters
    fitOptions.suppressFSPExpansion = true;
    DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
        DUSP1pars, fitOptions);
    ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars') 
end

%%    STEP 2.C. -- Plot predictions for other Dex concentrations.
showCases = [1,1,1,0];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1pars,Dusp1FitCases,showCases)
%%    STEP 2.D. -- Sample uncertainty for Dusp1 Parameters
%%      STEP 2.D.1. -- Compute sensitivity of the FSP solution
ModelGRDusp100nM.solutionScheme = 'fspSens';
sensSoln = ModelGRDusp100nM.solve();
ModelGRDusp100nM.solutionScheme = 'FSP';
%%      STEP 2.D.2. -- Compute FIM
% define which species in model are not observed.
ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};
% TODO - Make this automated when you load data.

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
%
%%      STEP 2.D.3. -- Run Metropolis Hastings Search
if loadPrevious
    MHDusp1File = 'MHDusp1_Jul22';
    load(MHDusp1File)
else
    MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
    MHFitOptions.thin=1;
    MHFitOptions.numberOfSamples=10000;
    MHFitOptions.burnIn=0;
    MHFitOptions.progress=true;
    MHFitOptions.numChains = 1;
    MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
    [DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
        [], MHFitOptions, 'MetropolisHastings');
    delete('TMPEricMHDusp1.mat')
    ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
end
%%      STEP 2.D.4. -- Plot the MH results
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1]; %
        CH(1).LineWidth = 3;
    end
end

%%
%%  STEP 3. -- Model Extensions using ODE Analyses
if loadPrevious
    vaNamesExtended = {'ModelGRfit'
        'extendedMod'
        'ModelGRDusp100nM_ext_red'
        'fullPars'
        };
    load(savedWorkspace,vaNamesExtended{:})

    try
        extendedMod.propensitiesGeneral{1}.stateDependentFactor(0)
    catch
        extendedMod = extendedMod.formPropensitiesGeneral('NucCyt_Model');
    end
    try
        ModelGRDusp100nM_ext_red.propensitiesGeneral{1}.stateDependentFactor(0)
    catch
        ModelGRDusp100nM_ext_red = ModelGRDusp100nM_ext_red.formPropensitiesGeneral('ExtModel100nm');
    end
    try
        for i=1:length(ModelGRfit)
            ModelGRfit{1}.propensitiesGeneral{1}.stateDependentFactor(0)
        end
    catch
        for i=1:length(ModelGRfit)
            ModelGRfit{i} = ModelGRfit{i}.formPropensitiesGeneral(['ExtModel100nm_',num2str(i)]);
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
    % However, we do not expec tthe FSP analysis to be very fast just yet, so
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
    extendedMod.parameters(4,:) = {'knuc2cyt',0.027};
    extendedMod.parameters(14,:) = {'knucdeg',0.001};
    extendedMod.parameters(15,:) = {'degCyt',0.01};
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
showCases = [1,1,1,0];
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
%%  STEP 4. -- Extend analysis to TS dynamics
% Load data that includes the TS.
allData = readtable('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv');
outlierThresh = 50;
fitGRinclude = true;

if loadPrevious
    vaNamesTS = {'ModelGRDusp100nM_ext_red'
        'parsAll_GR_Dusp1_TS'
        };
    load(savedWorkspace,vaNamesTS{:})
    if fitGRinclude
        % The set of parameters, including elongation time at the end.
        parsAllandTS = parsAll_GR_Dusp1_TS([1:12,14:16]);
    else
        parsAllandTS = parsAll_GR_Dusp1_TS([1:4,14:16]);
    end
else
    if fitGRinclude
        % The set of parameters, including elongation time at the end.
        parsAllandTS = [fullPars([1:12,13:14]),5];
    else
        parsAllandTS = [fullPars([1:4,13:14]),5];
    end
end

if fitGRinclude
    % The set of parameters, including elongation time at the end.
    indsDuspMod = [1:12,14];
    indsDuspDat = [1:13];
    indsTSmod = [1:12];
    indsTSpars = [1:12,15];
    indsODEmod = [1:12,14:15];
    indsODEdat = [1:14];
else
    indsDuspMod = [1:4,14];
    indsDuspDat = [1:4,5];
    indsTSmod = [1:3];
    indsTSpars = [1:3,7];
    indsODEmod = [1:4,14:15];
    indsODEdat = [1:6];
end
%%    STEP 4.A.1. -- Set up objective function to now include the TS analysis.
% Create prior for all parameters (just adding the elongation time to the
% previous list).
log10PriorMean = [-1 -1 0 -2,... %dusp1 pars
    -1 -3 -2 -1 -2 -2 -2 0.5, ...%GR pars
    NaN, ... % Dex concentration -- known
    -2, -3, ... % dusp1 transport, cyt RNA degradation
    0.7]; %elongation time (about 5 min)
log10PriorStd = 2*ones(1,16);

if fitGRinclude
    % Set parameters to be free during search.
    % This group allows the GR and DUSP1 parameters.
    ModelGRDusp100nM_ext_red.fittingOptions.modelVarsToFit = indsDuspMod;
    extendedMod.fittingOptions.modelVarsToFit = indsODEmod;
    
    for i = 1:3
        ModelGRfit{i}.tSpan = unique([0,ModelGRfit{i}.dataSet.times]);
    end

    Organization = {ModelGRfit{1},[5:12],[5:12],'computeLikelihood',1;...
        ModelGRfit{2},[5:12],[5:12],'computeLikelihood',1;...
        ModelGRfit{3},[5:12],[5:12],'computeLikelihood',1;...
        ModelGRDusp100nM_ext_red,indsDuspMod,indsDuspDat,'computeLikelihood',1;...
        extendedMod,indsODEmod,indsODEdat,'computeLikelihoodODE',0.01};
    logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:12,14,15,16])).^2./(2*log10PriorStd([1:12,14,15,16]).^2));
else
    % Set parameters to be free during search.
    % This group allows just the DUSP1 parameters to change and only fits
    % the DUSP1 data.
    ModelGRDusp100nM_ext_red.fittingOptions.modelVarsToFit = indsDuspMod;
    ModelGRDusp100nM_ext_red.fspOptions.fspTol = inf;
    extendedMod.fittingOptions.modelVarsToFit = indsODEmod;
    Organization = {ModelGRDusp100nM_ext_red,indsDuspMod,indsDuspDat,'computeLikelihood',1;...
        extendedMod,indsODEmod,indsODEdat,'computeLikelihoodODE',0.01};
    logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:4,14,15,16])).^2./(2*log10PriorStd([1:4,14,15,16]).^2));
end

% Create Copy of the prevous model
ModelTS = ModelGRDusp100nM_ext_red;
% This model can use all parameters except for the cytoplasmic degradation
ModelTS.fittingOptions.modelVarsToFit = indsTSmod;

% Objective function for just the TS analyses -- see function description
% below. 
objTS = @(x)computeTSlikelihood(x,ModelTS,allData,100,outlierThresh);

% Objective function for all analyses -- same as before but now with the
% transcription site analysis.
extraObjs(1,:) = {objTS,indsTSpars};
objAllwithTS = @(x)-getTotalFitErr(Organization,exp(x),false,extraObjs)-logPriorAll(exp(x));

%%    STEP 4.A.2. --  Test that the objective function works.
% This should take about 1s, and the value should be about 107700 with a
% good choice of parameters (when fitting all data including GR, DUSP1
% and TS). 
tic
    objAllwithTS(log(parsAllandTS))
toc
%%    STEP 4.B.1 -- Run the fit with fminsearch
fitOptions.MaxIter = 1000;
fitOptions.UseParallel = true;
for i = 1:8
    parsAllandTS = exp(fminsearch(objAllwithTS,log(parsAllandTS),fitOptions));
    bestObj = objAllwithTS(log(parsAllandTS));
end
%%    STEP 4.B.2 -- Run the fit again with particle swarm

% % This code could be used to run a particle swarm for optimization.  This
% works okay, but seems slower than the other approaches and did not reach
% as good an optimum as that found using the fminsearch routine.  I am
% leaving it here because it makes good use of parallel analyses, and might
% work better on a larger computer.
% LB = log(parsAllandTS)-log(20);
% UB = log(parsAllandTS)+log(20); 
% options = optimoptions('particleswarm', 'Display', 'iter');
% options.UseParallel = true;
% options.MaxIterations = 200; 
% options.SwarmSize = 200;
% X = exp(particleswarm(objAllwithTS,length(parsAllandTS),LB,UB,options));
% 
% Xobj = objAllwithTS(log(X))
% if Xobj<bestObj
%     parsAllandTS = X;
%     bestObj
% end


%%    STEP 4.C. -- Update the parameter set and make plots of results
%%      STEP 4.C.1. -- Make plots of the GR Dynamics.
if fitGRinclude
    makeGRPlots(combinedGRModel,parsAllandTS(5:12));
end
%%      STEP 4.C.2. -- Make plots of the Nucelar DUSP1 Dynamics.
ModelGRDusp100nM_ext_red.parameters(indsDuspMod,2) = num2cell(parsAllandTS(indsDuspDat));
showCases = [0,0,0,1];
makePlotsDUSP1({ModelGRDusp100nM_ext_red},ModelGRDusp100nM_ext_red,parsAllandTS(indsDuspDat),Dusp1FitCases,showCases)

%%      STEP 4.C.3. -- Make plots of the TS Dynamics.
TSthresh = 4;
ModelTS.parameters(indsTSmod,2) = num2cell(parsAllandTS(indsTSpars(1:end-1)));
% plot 100nm Fit Results
dex = 100;
[TSLikelihood,modelResults,dataResults] = computeTSlikelihood(parsAllandTS(indsTSpars),ModelTS,allData,dex,52,true); 
makeTSPlots(modelResults,dataResults,ModelTS,TSthresh,[700,701])

dex = 10;
[~,modelResults10,dataResults10] = computeTSlikelihood(parsAllandTS(indsTSpars),ModelTS,allData,dex,52,true);
dataResults10.meanNascentDat(1) = dataResults.meanNascentDat(1);
dataResults10.TS_counts(1) = dataResults.TS_counts(1);
dataResults10.fracTSDat(1) = dataResults.fracTSDat(1);
dataResults10.fracTSDatstd(1) = dataResults.fracTSDatstd(1);
dataResults10.meanNascentDatstd(1) = dataResults.meanNascentDatstd(1);
makeTSPlots(modelResults10,dataResults10,ModelTS,TSthresh,[702,703])

dex = 1;
[~,modelResults1,dataResults1] = computeTSlikelihood(parsAllandTS(indsTSpars),ModelTS,allData,dex,52,true); 
dataResults1.meanNascentDat(1) = dataResults.meanNascentDat(1);
dataResults1.TS_counts(1) = dataResults.TS_counts(1);
dataResults1.fracTSDat(1) = dataResults.fracTSDat(1);
dataResults1.fracTSDatstd(1) = dataResults.fracTSDatstd(1);
dataResults1.meanNascentDatstd(1) = dataResults.meanNascentDatstd(1);
makeTSPlots(modelResults1,dataResults1,ModelTS,TSthresh,[704,705])

for f = 701:2:705
    figure(f)
    subplot(2,1,1)
    set(gca,'ylim',[0,4])
    subplot(2,1,2)
    set(gca,'ylim',[0,0.5])
end

%%      STEP 4.C.4. -- Make plots of the TS Dynamics for titration
ModelTS75 = ModelTS;
ModelTS75.tSpan = [0,75];
dexV = unique(allData.Dex_Conc);
% modTitration
% datTitration
for i = 1:length(dexV)
    dex = dexV(i);
    AllDatRed = allData(allData.Time_index==75&allData.Dex_Conc==dex,:);
    [~,mod,dat] = computeTSlikelihood(parsAllandTS(indsTSpars),ModelTS75,AllDatRed,dex,52,true); 
    tmp.tSpan(i) = dex;
    modTitration.meanNascentMod(i) = mod.meanNascentMod(2);  
    modTitration.fracTSMod(i) = mod.fracTSMod(2);  
    modTitration.Psaved(i) = mod.Psaved(2); 

    datTitration.meanNascentDat(i) = dat.meanNascentDat(2);  
    datTitration.fracTSDat(i) = dat.fracTSDat(2);  
    datTitration.TS_counts(i) = dat.TS_counts(2);  
    
    % reps = unique(AllDatRed.Replica);
    % frac = [];
    % mn = [];
    % for j = 1:length(reps)
    %     AllDatRedReps = AllDatRed(strcmp(AllDatRed.Replica,reps{j}),:);
    %     frac(j) = sum(AllDatRedReps.DUSP1_ts_size_0>0)/size(AllDatRedReps,1);
    %     mn(j) = (sum(AllDatRedReps.DUSP1_ts_size_0(AllDatRedReps.DUSP1_ts_size_0>0))+...
    %          sum(AllDatRedReps.DUSP1_ts_size_1(AllDatRedReps.DUSP1_ts_size_1>0)))/size(AllDatRedReps,1);
    % end
    datTitration.fracTSDatstd(i) = dat.fracTSDatstd(2);
    datTitration.meanNascentDatstd(i) = dat.meanNascentDatstd(2);
end
makeTSPlots(modTitration,datTitration,tmp,TSthresh,[706,707])
figure(707)
subplot(2,1,1)
set(gca,'xscale','log','ylim',[0,4])
subplot(2,1,2)
set(gca,'xscale','log','ylim',[0,0.5])
xlabel('Dex Concentration (nM)')

%%      STEP 4.C.5. -- Make plots of the Cytoplasmic DUSP1 Dynamics.
extendedMod.parameters(indsODEmod,2) = num2cell(parsAllandTS(indsODEdat));

plotODEresults(extendedMod,extendedMod.solve,ModelGRfit{3},500)

extendedMod1 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
    {'rna','RNA_DUSP1_nuc'; ...
    'rCyt','RNA_DUSP1_cyto'},...
    {'Dex_Conc','1'});
extendedMod1.parameters(13,:) = {'Dex0',1.0};
plotODEresults(extendedMod1,extendedMod1.solve,ModelGRfit{1},501)

extendedMod10 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
    {'rna','RNA_DUSP1_nuc'; ...
    'rCyt','RNA_DUSP1_cyto'},...
    {'Dex_Conc','10'});
extendedMod10.parameters(13,:) = {'Dex0',10};
plotODEresults(extendedMod10,extendedMod10.solve,ModelGRfit{2},502)

%% Save Results for Easier Use in subsequent runs.
parsAll_GR_Dusp1_TS = [extendedMod.parameters{:,2}];
parsAll_GR_Dusp1_TS(indsODEmod) = parsAllandTS;
varNames = unique({'ModelGR',
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
    'parsAll_GR_Dusp1_TS'
    });

save('workspaceJuly22',varNames{:})

%%
%% Additonal Codes for FIM and PDO analyses (ALEX Paper)
%%  STEP XX -- Explore Optimal Designs Using FIM Add FIM to MH UQ Plots 
%%    STEP XX.A. -- FIM Analyses
%%      STEP XX.A.1. -- Plot UQ from FIM compared to MH
figNew = figure;

fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
    diag(log10PriorStd.^2));
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimTotal],'log',[],figNew)
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];
        CH(2).LineWidth = 3;
    end
end

%%      STEP XX.A.2. -- Find optimal experiment design (same number of cells)
nTotal = sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOpt = ModelGRDusp100nM.optimizeCellCounts(fimResults,nTotal,'tr[1:4]');
nCellsOptAvail = min(nCellsOpt,ModelGRDusp100nM.dataSet.nCells')
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(log10PriorStd.^2));
fimOptAvail = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptAvail,diag(log10PriorStd.^2));
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimTotal],'log',[],figNew);
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOptAvail,fimTotal],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];
        CH(2).LineWidth = 3;
        CH(3).Color=[0,1,1];
        CH(3).LineWidth = 3;
    end
end

f = figure;
set(f,'Position',[616   748   412   170])
bar([1:12],ModelGRDusp100nM.dataSet.nCells,0.45)
hold on
bar([1:12]+0.5,nCellsOpt',0.45)
set(gca,'xtick',[1:12]+0.25,'xticklabel',ModelGRDusp100nM.dataSet.times,'fontsize',16,'ylim',[0,7000])
legend('Intuitive Design','Optimal Design')


%%    STEP XX.B. -- PDO Calculations
%%      STEP XX.B.1. -- Calibrate PDO from Multi-Modal Experimental Data
% Calibration the PDO from empirical data. Here, the number of spots has
% been measured using different assays in data columns 'nTotal' for the
% 'true' data set and in the columns 'nSpots0' for a different label or
% 'intens1' for the integrated intensity.  We calibrate two different PDOs
% for this case. In both cases, we assume an 'AffinePoiss' PDO where the
% obervation probability is a Poisson distribution where the mean value is
% affine linearly related to the true value: P(y|x) = Poiss(a0 + a1*x);
ModelPDOSpots = ModelGRDusp100nM.calibratePDO('../../ExampleData/pdoCalibrationData.csv',...
    {'rna'},{'nTotal'},{'nSpots0'},'AffinePoiss',true);

%%      STEP XX.B.2. -- Calibrate PDO from Eric's DUSP1 Intensity Data.
ModelPDOIntensEric = ModelGRDusp100nM;
ModelPDOIntensEric = ModelPDOIntensEric.calibratePDO('EricData/pdoCalibrationData_EricIntensity_ConcHigh.csv',...
    {'rna'},{'RNA_DUSP1_nuc'},{'Nuc_DUSP1_avg_int_tot'},'AffinePoiss',true,[1,230,0.5]);

%%    STEP XX.C. -- FIM + PDO Analyses
%%      STEP XX.C.1. -- Analyze FIM with PDO for MCP/smFISH
fimsPDOSpot = ModelPDOSpots.computeFIM(sensSoln.sens,'log');
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(log10PriorStd.^2));

nCellsOptPDOspots = ModelPDOSpots.optimizeCellCounts(fimsPDOSpot,nTotal,'tr[1:4]');


figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimPDOSpots,fimTotal],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];
        CH(2).LineWidth = 3;
        CH(4).Color=[0,1,1];
        CH(4).LineWidth = 3;
        CH(3).Color=[0,0,1];
        CH(3).LineWidth = 3;
    end
end
%%      STEP XX.C.2. -- Analyze FIM with PDO for Intensity only
fimsPDOIntens = ModelPDOIntensEric.computeFIM(sensSoln.sens,'log');
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(log10PriorStd.^2));

nCellsOptPDOintens = ModelPDOSpots.optimizeCellCounts(fimsPDOIntens,nTotal,'tr[1:4]');

fimPDOIntensAvail = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptAvail,diag(log10PriorStd.^2));
figNew = figure; clf;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimPDOSpots,fimTotal,fimPDOIntens],'log',[],figNew);
% figNew = figure; clf;
% ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimTotal,fimPDOIntensAvail],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(3).Color=[0,0,0];
        CH(3).LineWidth = 3;
        CH(5).Color=[0,1,1];
        CH(5).LineWidth = 3;
        CH(4).Color=[0,0,1];
        CH(4).LineWidth = 3;
        CH(2).Color=[0,1,0];
        CH(2).LineWidth = 3;
    end
end

%%      STEP XX.C.3. -- Analyze FIM with PDO for Intensity only more cells
nTimes = 3.71;
fimPDOIntens2x = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt*nTimes,diag(log10PriorStd.^2));
det(fimOpt{1}(1:4,1:4))/det(fimPDOIntens2x{1}(1:4,1:4))
figNew = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimPDOSpots,fimTotal,fimPDOIntens,fimPDOIntens2x],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(4).Color=[0,0,0];
        CH(4).LineWidth = 3;
        CH(6).Color=[0,1,1];
        CH(6).LineWidth = 3;
        CH(5).Color=[0,0,1];
        CH(5).LineWidth = 3;
        CH(3).Color=[0,1,0];
        CH(3).LineWidth = 3;
        CH(2).Color=[1,0,0];
        CH(2).LineWidth = 3;
    end
end



%%    STEP XX.D. -- Validation of FIM Predictions
%%      STEP XX.D.1. -- Fit to FIM selected time points
ModelGRDusp100nM_FIMDesign = ModelGRDusp100nM;
% Set the fitting routine to only consider the time points selected by the
% FIM analysis:
ModelGRDusp100nM_FIMDesign.fittingOptions.timesToFit = nCellsOpt>0;
% Refit the model, but now with only those time points.
% DUSP1parsFIMDesign = DUSP1pars;
for i = 1:1
    DUSP1parsFIMDesign = ModelGRDusp100nM_FIMDesign.maximizeLikelihood(...
        DUSP1parsFIMDesign, fitOptions);
    ModelGRDusp100nM_FIMDesign.parameters(1:4,2) = num2cell(DUSP1parsFIMDesign);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign') 
end

%%      STEP XX.D.2. -- Plot fit and predicitons using FIM suggested conditions.
showCases = [1,1,1,0];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsFIMDesign,Dusp1FitCases,showCases)

%%      STEP XX.D.3. -- Fit the DISTORTED intensity measurements at ALL times.
% Here, I will load the intensity data for the nuclear DUSP1.  then I will
% attempt to identify the model from just that data and at just the times
% selected by the FIM.
ModelPDOIntensEric.dataSet = [];
ModelPDOIntensEric = ModelPDOIntensEric.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','Nuc_DUSP1_avg_int_tot'},...
        {'Dex_Conc','100'}); 
load('EricModelDusp1_MMDex','DUSP1parsIntensity') 

ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1parsIntensity);

for i = 1:fitIters
    DUSP1parsIntensity = ModelPDOIntensEric.maximizeLikelihood(...
        DUSP1parsIntensity, fitOptions);
    ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1parsIntensity);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsIntensity','DUSP1parsFIMDesignIntensity') 
end
%%        STEP XX.D.3.a. -- Plot predictions when fit to distorted data at ALL times.
showCases = [1,1,1,0];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsIntensity,Dusp1FitCases,showCases)

%%      STEP XX.D.4. -- Fit the DISTORTED intensity measurements at FIM selected times.
% Here, I will load the intensity data for the nuclear DUSP1.  then I will
% attempt to identify the model from just that data and at just the times
% selected by the FIM.
ModelPDOIntensEricFIM = ModelPDOIntensEric;
ModelPDOIntensEricFIM.dataSet = [];
ModelPDOIntensEricFIM = ModelPDOIntensEricFIM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','Nuc_DUSP1_avg_int_tot'},...
        {'Dex_Conc','100'}); 
% Set fitting routine only to consider the time points selected by the FIM.
ModelPDOIntensEricFIM.fittingOptions.timesToFit = nCellsOpt>0;
% Refit the model, but now with only those time points.
load('EricModelDusp1_MMDex','DUSP1parsFIMDesignIntensity') 
ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesignIntensity);

for i = 1:5
    DUSP1parsFIMDesignIntensity = ModelPDOIntensEricFIM.maximizeLikelihood(...
        DUSP1parsFIMDesignIntensity, fitOptions);
    ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesignIntensity);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsFIMDesignIntensity') 
end
%%        STEP XX.D.5.a. -- Plot predictions when fit to distorted data at FIM times.
showCases = [1,1,1,0];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsFIMDesignIntensity,Dusp1FitCases,showCases)

%%    STEP XX.E. -- Plot Information vs. Expt Design, PDO, and Number of Cells
fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,diag(log10PriorStd.^2));
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(log10PriorStd.^2));
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(log10PriorStd.^2));
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(log10PriorStd.^2));

barWithOriginalNumber = [det(fimTotal{1}(1:4,1:4)),det(fimOpt{1}(1:4,1:4)),det(fimPDOSpots{1}(1:4,1:4)),det(fimPDOIntens{1}(1:4,1:4))];

nCellVec = logspace(3,5,20);
nCellsOrigRat = ModelGRDusp100nM.dataSet.nCells/sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOptRat = nCellsOpt/sum(nCellsOpt);

cols = {'k','c','b','r'};
fimDetVsNumberMAt=[];
for i = 1:length(nCellVec)
    fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOrigRat*nCellVec(i),diag(log10PriorStd.^2));
    fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptRat*nCellVec(i),diag(log10PriorStd.^2));
    fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOptRat*nCellVec(i),diag(log10PriorStd.^2));
    fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptRat*nCellVec(i),diag(log10PriorStd.^2));
    fimDetVsNumberMAt(i,:) = [det(fimOrig{1}(1:4,1:4)),det(fimOpt{1}(1:4,1:4)),det(fimPDOSpots{1}(1:4,1:4)),det(fimPDOIntens{1}(1:4,1:4))];
end

close all
figNew = figure;
% loglog(nCellVec,1./fimDetVsNumberMAt,'linewidth',2);
for i = 1:4
    loglog([10,1e6],1/barWithOriginalNumber(i)*[1,1],[cols{i},'--'],'linewidth',2)
    hold on
    loglog(nCellVec,1./fimDetVsNumberMAt(:,i),cols{i},'linewidth',2)
end
set(gca,'fontsize',15,'ylim',10.^[-12.,-9],'xlim',10.^[3,5])


%% Extra Functions
function makeGRPlots(combinedModel,GRpars)
combinedGRModel = combinedModel.updateModels(GRpars,false);
nMods = length(combinedGRModel.SSITModels);
ModelGroup = cell(nMods,1);
for i=1:nMods
    %  Update parameters in original models.
    ModelGroup{i} = combinedGRModel.SSITModels{i};
    ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
    ModelGroup{i}.makeFitPlot([],1,[],true,'STD')
end
end

function plotODEresults(extendedMod,soln,modeWithGRData,fignum)
arguments
    extendedMod
    soln
    modeWithGRData
    fignum = 1;
end
figure(fignum); clf;
% Plot GR levels vs. Time
subplot(2,1,1)
plot(extendedMod.tSpan,soln.ode(:,3:4),'--','LineWidth',2);hold on
plot(modeWithGRData.dataSet.times,modeWithGRData.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
legend(extendedMod.species(3:4))
set(gca,'xlim',[-10,200],'ylim',[0,12],'fontsize',16)
ylabel('GR Concentrations (UA)')
legend({'Cyt-Model','Nuc-Model','Cyt-Data','Nuc-Data'})
title('GR')

% Plot DUSP1 levels vs. Time
subplot(2,1,2)
plot(extendedMod.tSpan,soln.ode(:,5:6),'--','LineWidth',2);hold on
plot(extendedMod.dataSet.times,extendedMod.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
legend(extendedMod.species(3:4))
set(gca,'xlim',[-10,200],'ylim',[0,160],'fontsize',16)
ylabel('GR Concentrations (UA)')
legend({'Nuc-Model','Cyt-Model','Nuc-Data','Cyt-Data'})
title('DUSP1')
end

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

function [TSLikelihood,modelResults,dataResults] = computeTSlikelihood(x,ModelTS,allData,dex,outlierThresh,reportErrors)
arguments
    x
    ModelTS
    allData
    dex = 100
    outlierThresh = 50;  
    reportErrors = false
    % The outlierThresh input sets an outlier threshold for the nascent RNA
    % at the TS.  Data with observed counts more than this value will be
    % trincated to this value.  Model predictions for counts higher than
    % this value will also be integrated into a single bin "more than
    % outlierThresh".
end

tau_elong = x(end);  % Time for elongation and processing of mRNA (assume deterministic)

% The rest of the parameters go into the model
x = x(1:end-1);
ModelTS.parameters(ModelTS.fittingOptions.modelVarsToFit,2) = num2cell(x);
ModelTS.parameters(13,:) = {'Dex0',dex};

% set upper bound on FSP solution
ModelTS.fspOptions.bounds(6) = outlierThresh*1.5;

% turn off FSP expansion routine
ModelTS.fspOptions.fspTol = inf;

% Create past model only used to get distribution of states for when
% elongation would begin.
ModelTSPast = ModelTS;
ModelTSPast.parameters(3,:) = {'kr',0};  % Turn off transcription.
ModelTSPast.initialCondition(5) = 0;

% Shift solution time points back by tau_elong.
ModelTSPast.tSpan = ModelTS.tSpan - tau_elong;
ModelTSPast.tSpan = [ModelTSPast.tSpan(1)-10,ModelTSPast.tSpan];
ModelTSPast.initialTime = ModelTSPast.tSpan(1);

% Solve for distributions at t - tau_elong.
TSPastSoln = ModelTSPast.solve;

% For the prsent model we start one elongation period in the past and then 
% integrate forward in time.
ModelTSPresent = ModelTS;
ModelTSPresent.parameters(4,:) = {'knuc2cyt',0};
ModelTSPresent.parameters(14,:) = {'knucdeg',0};

% We will NOT use steady state initial conditions -- rather, we use the
% computed distributions from previous time step.
ModelTSPresent.fspOptions.initApproxSS = false;
% We need to solve once for each time point, where the initial condition is
% defined by the past model.
NT = length(TSPastSoln.fsp)-1;
fsp = cell(1,NT);
for iT = 1:NT
    tmpModelTSPresent = ModelTSPresent.setICfromFspVector(TSPastSoln.stateSpace,TSPastSoln.fsp{iT+1});
    tmpModelTSPresent.tSpan = ModelTSPast.tSpan(iT+1)+[0,tau_elong/2,tau_elong];
    tmpModelTSPresent.initialTime = tmpModelTSPresent.tSpan(1);
    tsPresentSoln = tmpModelTSPresent.solve;
    fsp(iT) = tsPresentSoln.fsp(end); 
end

% Summarize the model results for plotting and for likelihood function
% calculations.
modelResults.meanNascentMod = zeros(NT,1);
modelResults.fracTSMod = zeros(NT,1);
TSthresh = 4;
for iT = NT:-1:1
    P=double(fsp{iT}.p.sumOver([1,2]).data);
    P = max(0,P);
    if sum(P)>1
        P = P/sum(P);
    end
    P(1:TSthresh) = [sum(P(1:TSthresh));zeros(TSthresh-1,1)]; 
    modelResults.Psaved{iT} = P;
end

% Summarize the data for plotting and for likelihood function
% calculations.
for iT = NT:-1:1
    time = ModelTS.tSpan(iT);
    redData = allData(allData.Dex_Conc==dex&allData.Time_index==time,:);
    TS_counts0 = redData.DUSP1_ts_size_0;
    TS_counts1 = redData.DUSP1_ts_size_1;    
    TS_counts0(isnan(TS_counts0)) = 0;
    TS_counts1(isnan(TS_counts1)) = 0;
    dataResults.TS_counts{iT} = TS_counts0+TS_counts1; 
    dataResults.meanNascentDat(iT) = mean(dataResults.TS_counts{iT});
    dataResults.fracTSDat(iT) = sum(dataResults.TS_counts{iT}>=TSthresh)/length(dataResults.TS_counts{iT});

    if reportErrors
        reps = unique(redData.Replica);
        frac = [];
        mn = [];
        for j = 1:length(reps)
            AllDatRedReps = redData(strcmp(redData.Replica,reps{j}),:);
            frac(j) = sum(AllDatRedReps.DUSP1_ts_size_0>0)/size(AllDatRedReps,1);
            mn(j) = (sum(AllDatRedReps.DUSP1_ts_size_0(AllDatRedReps.DUSP1_ts_size_0>0))+...
                sum(AllDatRedReps.DUSP1_ts_size_1(AllDatRedReps.DUSP1_ts_size_1>0)))/size(AllDatRedReps,1);
        end
        dataResults.fracTSDatstd(iT) = std(frac);
        dataResults.meanNascentDatstd(iT) = std(mn);
    end
end


% Adjust to account for partial transcripts.
nRNAmax = floor(ModelTS.fspOptions.bounds(6)+1);
PDO_partial_transcripts = zeros(nRNAmax,nRNAmax);
for i = 1:nRNAmax
    PDO_partial_transcripts(1:i,i) = binopdf(0:i-1,i-1,0.5);
end

% Calculate the likelihood function from the TS data.
TSLikelihood=0;
for iT=1:NT
    modelResults.Psaved{iT} = max(1e-10,modelResults.Psaved{iT});
    
    % Remove outliers from data and set to threshold.
    countsRNA = dataResults.TS_counts{iT};
    countsRNA(countsRNA>outlierThresh) = outlierThresh;
   
    % Apply PDO to account for partial transcripts.
    try
        modelResults.Psaved{iT} = PDO_partial_transcripts*modelResults.Psaved{iT};
    catch
        1+1
    end
    % Bin outlier predictions for model
    P = modelResults.Psaved{iT};

    P(outlierThresh+1) = sum(P(outlierThresh+1:end));
    P = P(1:outlierThresh+1);

    modelResults.meanNascentMod(iT) = [TSthresh:length(P)-1]*P(TSthresh+1:end);
    modelResults.fracTSMod(iT) = sum(P(TSthresh+1:end));

    TSLikelihood = TSLikelihood + sum(log(P(countsRNA+1)));
end
end

function makeTSPlots(modelResults,dataResults,ModelTS,TSthresh,figNum)
arguments
    modelResults
    dataResults
    ModelTS
    TSthresh = 4;
    figNum = [700,701];
end
figure(figNum(1)); clf
NT = length(modelResults.meanNascentMod);
for iT = 1:NT
    P=modelResults.Psaved{iT};
   
    % Lump predictions less than detection limit at zero.
    P(1:TSthresh) = sum(P(1:TSthresh))/4;  
    % P(1:TSthresh) = [sum(P(1:TSthresh));zeros(TSthresh-1,1)];   
    
    subplot(4,3,iT)
    % stairs([0:length(P)-1],cumsum(P),'linewidth',3)
    stairs([0:length(P)-1],(P),'linewidth',3)
    set(gca,'fontsize',16,'ylim',[0,1],'xlim',[0,50])
    set(gca,'fontsize',16,'ylim',[0,0.02],'xlim',[0,50])
    if iT>9
        xlabel('Number nascent RNA')
    end
    if mod(iT-1,3)==0
        ylabel('Probability')
    end   
end

figure(figNum(2)); clf
subplot(2,1,1)
plot(ModelTS.tSpan,modelResults.meanNascentMod,'linewidth',3);
set(gca,'fontsize',16,'ylim',[0,3])
ylabel('Mean Nascent RNA')

subplot(2,1,2)
plot(ModelTS.tSpan,modelResults.fracTSMod,'linewidth',3);
set(gca,'fontsize',16,'ylim',[0,0.3])
ylabel('Fraction with TS')
xlabel('time (min)')
%
%    STEP 4.C.2. -- Add data to plots
figure(figNum(1))
for iT = 1:NT
    subplot(4,3,iT); hold on
    histogram(dataResults.TS_counts{iT},[0,TSthresh:52],'Normalization','pdf','DisplayStyle','bar','LineWidth',1);
    % histogram(dataResults.TS_counts{iT},[0,TSthresh:52],'Normalization','cdf','DisplayStyle','stairs','LineWidth',3);
    % nCount(iT) = length(dataResults.TS_counts{iT});
    % semCount(iT) = std(dataResults.TS_counts{iT})/sqrt(nCount(iT)-1);
end

figure(figNum(2))
subplot(2,1,1)
hold on
% plot(ModelTS.tSpan,dataResults.meanNascentDat,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
errorbar(ModelTS.tSpan,dataResults.meanNascentDat,dataResults.meanNascentDatstd,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');

subplot(2,1,2)
hold on
% stdFrac = sqrt(nCount.*dataResults.fracTSDat.*(1-dataResults.fracTSDat))./nCount;
% plot(ModelTS.tSpan,dataResults.fracTSDat,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
errorbar(ModelTS.tSpan,dataResults.fracTSDat,dataResults.fracTSDatstd,'s','markersize',16,'linewidth',3,'MarkerFaceColor','k');
end