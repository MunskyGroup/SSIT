%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations
close all 
addpath(genpath('../../src'));

GR_Data = 'RonData062025/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv';
Dusp1Data = 'RonData062025/ssit_gated_Jun19';

loadPrevious = false;
savedWorkspace = 'workspaceJune25';
addpath('tmpPropensityFunctions');

fitOptions = optimset('Display','iter','MaxIter',10);
fitIters = 1;
makePlots = false;

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

    try
        ModelGRfit{1}.propensitiesGeneral{1}.stateDependentFactor(0);
    catch
        for i=1:length(ModelGRfit)
            ModelGRfit{i} = ModelGRfit{i}.formPropensitiesGeneral('GR_Model');
        end
    end
else
    %% STEP 0.B.1. -- Create Base Model for GR Only
    % Here, I set up the model for the GR translocation dynamics.
    ModelGR = SSIT;
    ModelGR.species = {'cytGR';'nucGR'};
    ModelGR.initialCondition = [20;1];

    ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';...
        'kg1';'gg1*cytGR';'gg2*nucGR'};
    ModelGR.stoichiometry = [-1,1,1,-1,0;...
        1,-1,0,0,-1];

    ModelGR.parameters = {'koff',0.1;'kon',0.1;'kr',1;'knuc2cyt',0.02;...
        'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',14e-5;...
        'gg1',1e-5;'gg2',1e-6;'MDex',5;'Dex0',100;...
        'knucdeg',0.0001;'degCyt',0.05;'degCytA',0.01};%'ktranslation',1;'gprotein',0.003;'degCyt2',1e-5});

    log10PriorMean = [-1 -1 0 -2,...
        -1 -3 -2 -1 -2,...
        -2 -2 0.5, NaN, ...
        -4, -1, 2];
    log10PriorStd = 2*ones(1,16);
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
        % ModelGRfit{i} = ModelGR.loadData(GR_Data,...
        %     {'nucGR','normGRnuc';'cytGR','normGRcyt'},...
        %     {'dex_conc',GRfitCases{i,2}});
        if i==3
            ModelGRfit{i} = ModelGR.loadData(GR_Data,...
                {'cytGR','normGRcyt';'nucGR','normGRnuc'},...
                {[],[], ...
                ['(TAB.dex_conc==',GRfitCases{i,1},'|TAB.dex_conc==0)&TAB.time~=20&TAB.time~=40&TAB.time~=60&TAB.time~=90&TAB.time~=150']});
        else
            ModelGRfit{i} = ModelGR.loadData(GR_Data,...
                {'cytGR','normGRcyt';'nucGR','normGRnuc'},...
                {[],[], ...
                ['TAB.dex_conc==',GRfitCases{i,1},'&TAB.time~=20&TAB.time~=40&TAB.time~=60&TAB.time~=90&TAB.time~=150']});
        end

        ModelGRfit{i}.parameters(13,:) = {'Dex0', str2num(GRfitCases{i,1})};
        ModelGRparameterMap(i) = {(1:8)};
        % parameters 1 - 8 refer to the parameter set that is relevant to
        % the entire class of models.  In this case, these refer to
        % global parameters 5:15 (GR parameters).
    end
    %% STEP 0.B.4. -- Make Guesses for the FSP bounds
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
%%
%% STEP 1 -- GR Model
%%     STEP 1.A. -- Combine all three GR models and fit using a single parameter set.
for i = 1:3
    ModelGRfit{i}.tSpan = ModelGRfit{i}.dataSet.times;
end
 
logPriorGR = @(x)-sum((log10(x)-log10PriorMean(5:12)).^2./(2*log10PriorStd(5:12).^2),"all");
for jj = 1:fitIters
    combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap,logPriorGR);
    combinedGRModel = combinedGRModel.initializeStateSpaces(boundGuesses);
    combinedGRModel = combinedGRModel.updateModels(GRpars,false);
    GRpars = combinedGRModel.maximizeLikelihood( GRpars, fitOptions);
    save('EricModel_MMDex','GRpars') 
end
% Best error so far: 33804.7
%%     STEP 1.B. -- Compute FIM for GR parameters.
% combinedGRModel = combinedGRModel.computeFIMs([],'log');
% fimGR_withPrior = combinedGRModel.FIM.totalFIM+... % the FIM in log space.
%     diag(1./(log10PriorStd(ModelGR.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.

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
if makePlots
    figNums = [1:12];
    for i = figNums
        try
            close(i);
        catch
        end
    end
    makeGRPlots(combinedGRModel,GRpars,GR_Data,false)
end
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
    ModelGRDusp100nM.customConstraintFuns = [];
    ModelGRDusp100nM.species = {'offGene';'onGene';'cytGR';'nucGR';'rna'};
    ModelGRDusp100nM.initialCondition = [2;0;24;1;5];
    ModelGRDusp100nM.propensityFunctions = {'(kon*nucGR)*offGene';'koff*onGene';
        '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR';...
        'kr*onGene';'knuc2cyt*rna'};
    ModelGRDusp100nM.stoichiometry = [-1,1,0,0,0,0,0,0,0;...
        1,-1,0,0,0,0,0,0,0;...
        0,0,-1,1,1,-1,0,0,0;...
        0,0,1,-1,0,0,-1,0,0;...
        0,0,0,0,0,0,0,1,-1];
    ModelGRDusp100nM.useHybrid = true;
    ModelGRDusp100nM.hybridOptions.upstreamODEs = {'cytGR','nucGR'};
    ModelGRDusp100nM.solutionScheme = 'FSP';
    ModelGRDusp100nM.fspOptions.bounds = [0;0;0;2;2;500];
    ModelGRDusp100nM.fittingOptions.modelVarsToFit = [1:4];
    ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('EricModDusp1');
    duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(ModelGRDusp100nM.fittingOptions.modelVarsToFit)).^2./(2*log10PriorStd(ModelGRDusp100nM.fittingOptions.modelVarsToFit).^2));
    ModelGRDusp100nM.fittingOptions.logPrior = duspLogPrior;

    %% STEP 2.A.2. -- Load pre-fit parameters into model.
    load('EricModelDusp1_MMDex','DUSP1pars')
    ModelGRDusp100nM.parameters(ModelGRDusp100nM.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);

    %% STEP 2.A.3. -- Load and Associate with DUSP1 smFISH Data (100nM Dex Only)
    % The commented code below would be needed to fit multiple conditions,
    % but that is not used in this case.  It is left here in case it is
    % needed in later stages of the project.
    Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
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

    % ModelGRDusp100nM = ModelGRDusp100nM.loadData(Dusp1Data,...
    %     {'rna','totalNucRNA'},{'Dex_Conc','100'});
    ModelGRDusp100nM = ModelGRDusp100nM.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'},...
        {[],[],['(TAB.dex_conc==100|TAB.dex_conc==0)', ...
        '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});

    DUSP1pars = [ModelGRDusp100nM.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
end

%%    STEP 2.B. -- Fit DUSP1 model at 100nM Dex.
for i = 1:fitIters
    fitOptions.suppressFSPExpansion = true;
    DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
        DUSP1pars, fitOptions);
    ModelGRDusp100nM.parameters(ModelGRDusp100nM.fittingOptions.modelVarsToFit,2) = ...
        num2cell(DUSP1pars);
    save('EricModelDusp1_MMDex','GRpars','DUSP1pars') 
end
% Best error so far: 17974.36

%%    STEP 2.C. -- Plot predictions for other Dex concentrations.
if makePlots
    figNums = [201,301,302,303,221,321,322,323];
    for i = figNums
        try
            close(i);
        catch
        end
    end
    showCases = [1,1,1,0];
    makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1pars,Dusp1FitCases,showCases,Dusp1Data);
end
%%    STEP 2.D. -- Sample uncertainty for Dusp1 Parameters
%%      STEP 2.D.1. -- Compute sensitivity of the FSP solution
% ModelGRDusp100nM.solutionScheme = 'fspSens';
% sensSoln = ModelGRDusp100nM.solve();
% ModelGRDusp100nM.solutionScheme = 'FSP';
%%      STEP 2.D.2. -- Compute FIM
% % define which species in model are not observed.
% ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};
% % TODO - Make this automated when you load data.
% 
% % compute the FIM
% fimResults = ModelGRDusp100nM.computeFIM(sensSoln.sens,'log');
% 
% % In the following, the log-prior is used as a prior co-variance matrix.
% % This will be used in the FIM calculation as an FIM without new evidence 
% % being set equal to the inverse of this covariance matrix.  More rigorous
% % justification is needed to support this heuristic.
% fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
%     diag(log10PriorStd(1:19).^2));
% 
% FIMfree = fimTotal{1}(1:4,1:4);
% if min(eig(FIMfree))<1
%     disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
%     FIMfree = FIMfree + 1*eye(length(FIMfree));
% end
% covFree = FIMfree^-1;
% covFree = 0.5*(covFree+covFree');
%
%%      STEP 2.D.3. -- Run Metropolis Hastings Search
% if loadPrevious
%     MHDusp1File = 'MHDusp1_Jul22';
%     load(MHDusp1File)
% else
%     MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
%     MHFitOptions.thin=1;
%     MHFitOptions.numberOfSamples=10000;
%     MHFitOptions.burnIn=0;
%     MHFitOptions.progress=true;
%     MHFitOptions.numChains = 1;
%     MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
%     [DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
%         [], MHFitOptions, 'MetropolisHastings');
%     delete('TMPEricMHDusp1.mat')
%     ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
% end
% %%      STEP 2.D.4. -- Plot the MH results
% figNew = figure;
% ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
% for i = 1:3
%     for j = i:3
%         subplot(3,3,(i-1)*3+j)
%         CH = get(gca,'Children');
%         CH(1).Color=[1,0,1]; %
%         CH(1).LineWidth = 3;
%     end
% end

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
    % and protein
    extendedMod = ModelGRDusp100nM;
    extendedMod = extendedMod.addSpecies({'rCyt'},0);
    % extendedMod = extendedMod.addSpecies({'pCyt'},0);
    % extendedMod = extendedMod.addSpecies({'pJNK'},0);

    % Adjust the reaction to be nuc -> cyt transport.
    % extendedMod.propensityFunctions{9,1} = 'knuc2cyt*rna';
    % extendedMod.parameters{4,1} = 'knuc2cyt';
    extendedMod.stoichiometry(:,9) = 0;
    extendedMod.stoichiometry([5,6],9) = [-1;1];

    % Add nuc degradation reaction.
    % extendedMod.propensityFunctions{10,1} = 'knucdeg*rna';
    % extendedMod.stoichiometry(:,10) = 0;
    % extendedMod.stoichiometry(5,10) = -1;
    % extendedMod.parameters(14,:) = {'knucdeg',1.4e-5};

    % Add cytoplasmic degradation as a reaction
    % extendedMod.propensityFunctions{10,1} = 'degCyt*rCyt';
    extendedMod.stoichiometry(:,10) = 0;
    extendedMod.stoichiometry(6,10) = -1;
    % extendedMod.parameters(15,:) = {'degCyt',0.036};
    % extendedMod.parameters(16,:) = {'degCytA',1/500};
    extendedMod.propensityFunctions{10,1} = 'degCyt.*(1+1./(1+degCytA*rCyt)).*rCyt';

    % % Add protein translation as a reaction
    % extendedMod.propensityFunctions{12,1} = 'ktranslation*rCyt';
    % extendedMod.stoichiometry(:,12) = 0;
    % extendedMod.stoichiometry(7,12) = 1;
    % extendedMod.parameters(16,:) = {'ktranslation',1e-9};

    % % Add protein degradation as a reaction
    % extendedMod.propensityFunctions{13,1} = 'gprotein*pCyt';
    % extendedMod.stoichiometry(:,13) = 0;
    % extendedMod.stoichiometry(7,13) = -1;
    % extendedMod.parameters(17,:) = {'gprotein',1};

    % % Add JNK phosphorylation
    % extendedMod.propensityFunctions{14,1} = 'kJNK';
    % extendedMod.stoichiometry(:,14) = 0;
    % extendedMod.stoichiometry(8,14) = 1;
    % extendedMod.parameters(18,:) = {'kJNK',1};
    % 
    % % Add JNK de-phosphorylation
    % extendedMod.propensityFunctions{15,1} = 'gJNK*pJNK.*pCyt';
    % extendedMod.stoichiometry(:,15) = 0;
    % extendedMod.stoichiometry(8,15) = -1;
    % extendedMod.parameters(19,:) = {'gJNK',1};
    % 
    % % Adjust other reactions based on protein concentration.
    % extendedMod.parameters(20,:) = {'fJNK',1e-5};
    % extendedMod.propensityFunctions{2,1} = 'koff./(1+fJNK*pJNK).*onGene';
    % % change to increase gene deactivation

    log10PriorMean = [-1 -1 0 -2,...
        -1 -3 -2 -1 -2,...
        -2 -2 0.5, NaN, ...
        -4, -1, 0,-3,-5,0,-5];
    log10PriorStd = 2*ones(1,20);

    % extendedMod.propensityFunctions{1,1} = '(kon*nucGR)/(1+pCyt)*offGene';
    % change to decrease gene activation

    % % The following code can be used to verify that the model changes so far
    % % are consistent with the previous FSP model.
    % % This is useful to verify before changing the reactions too much.
    % % However, we do not expec tthe FSP analysis to be very fast just yet, so
    % % it would not make much sense to run this for fitting.
    % if false
    %     extendedMod.initialCondition(6) = 1;
    %     extendedMod.fspOptions.bounds = [0;0;0;0;2;2;400;1]
    %     extendedMod.fspOptions.verbose = 1;
    %     extendedMod = extendedMod.formPropensitiesGeneral('extendedModel');
    %     [~,extendedMod.fspOptions.bounds] = extendedMod.solve;
    %     extendedMod = extendedMod.loadData(Dusp1Data,...
    %         {'rna','num_nuc_spots'; ...
    %         'rCyt','num_cyto_spots'},...
    %         {'dex_conc','100'});
    %     extendedMod.makeFitPlot
    % end
    %%
    
    %%
    % ModelCytOnly = extendedMod;
    % ModelCytOnly.hybridOptions.upstreamODEs = {'offGene','onGene','cytGR','nucGR','rna'};
    % ModelCytOnly = ModelCytOnly.formPropensitiesGeneral('ODEEricCytOnly');
    % % ModelCytOnly.customConstraintFuns = {'abs(rna-rCyt)'};
    % ModelCytOnly.fspOptions.verbose = true;
    % % ModelCytOnly.fspOptions.bounds = [0,0,0,0,2,2,288,300,100];
    % ModelCytOnly.fspOptions.bounds = [0,300];
    % ModelCytOnly.parameters{15,2} = 0.01;
    % [solnCytOnly,~,ModelCytOnly] = ModelCytOnly.solve;
    % 
    % %   STEP 3.B.1 -- Load data into the model
    % dex = 100; dexModel = 3;
    % ModelCytOnly.parameters(13,:) = {'Dex0',dex};
    % ModelCytOnly = ModelCytOnly.loadData(Dusp1Data,...
    %     {'rCyt','num_cyto_spots'},...
    %     {[],[],['(TAB.dex_conc==',num2str(dex),'|TAB.dex_conc==0)', ...
    %     '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});

    %%
    % ModelCytOnly.fittingOptions.modelVarsToFit = [4,14:16];
    % ModelCytOnly.fittingOptions.logPrior = ...
    %     @(x)-sum((log10(x(:))'-log10PriorMean(ModelCytOnly.fittingOptions.modelVarsToFit)).^2./...
    %     (2*log10PriorStd(ModelCytOnly.fittingOptions.modelVarsToFit).^2));
    % 
    % pars = [ModelCytOnly.parameters{ModelCytOnly.fittingOptions.modelVarsToFit,2}];
    % for i = 1:fitIters
    %     pars = ModelCytOnly.maximizeLikelihood(pars,fitOptions);
    % end
    % ModelCytOnly.parameters(ModelCytOnly.fittingOptions.modelVarsToFit,2) = num2cell(pars);
    % ModelCytOnly.makeFitPlot

    %%    STEP 3.A.2 -- Create reduced model for Nuclear DUSP1 but with same parameters as new model
    % This model does not have the cytoplasmic RNA, protein or feedback, so
    % it is expected only to be valid for early times.
    ModelGRDusp100nM_ext_red = ModelGRDusp100nM;

    % Adjust the final reaction to be nuc -> cyt transport.
    ModelGRDusp100nM_ext_red.propensityFunctions{9,1} = '(knuc2cyt+knucdeg)*rna';
    % ModelGRDusp100nM_ext_red.parameters{4,1} = 'knuc2cyt';
    % ModelGRDusp100nM_ext_red.stoichiometry(:,9) = 0;
    % ModelGRDusp100nM_ext_red.stoichiometry(5,9) = -1;

    % Add nuc degradation reaction.
    % ModelGRDusp100nM_ext_red.propensityFunctions{10,1} = 'knucdeg*rna';
    % ModelGRDusp100nM_ext_red.stoichiometry(:,10) = 0;
    % ModelGRDusp100nM_ext_red.stoichiometry(5,10) = -1;
    % ModelGRDusp100nM_ext_red.parameters(14,:) = {'knucdeg',1.5e-5};

    % This model does not have the cytoplasmic RNA, protein or feedback, so
    % it is expected only to be valid for early times. We will remove all
    % but the first few time points for greater efficiency.
    ModelGRDusp100nM_ext_red.fittingOptions.timesToFit = [1:2];
    ModelGRDusp100nM_ext_red.tSpan = ModelGRDusp100nM_ext_red.tSpan(1:2);

    ModelGRDusp100nM_ext_red = ModelGRDusp100nM_ext_red.formPropensitiesGeneral('ModelGRDusp100nM_ext_red');
    ModelGRDusp100nM_ext_red.fittingOptions.modelVarsToFit = [1:4,14];
    ModelGRDusp100nM_ext_red.fittingOptions.logPrior = [];
    if makePlots
        ModelGRDusp100nM_ext_red.makeFitPlot
    end
    % This just makes another plot to verify that the model is still
    % matching for the extra reaction.

    %%    STEP 3.B.1 -- Load data into the model
    dex = 100; dexModel = 3;
    extendedMod.parameters(13,:) = {'Dex0',dex};
    extendedMod = extendedMod.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],['(TAB.dex_conc==',num2str(dex),'|TAB.dex_conc==0)', ...
        '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});

    %%    STEP 3.B.2 -- Switch solver to ODE and generate model codes
    extendedMod.solutionScheme = 'ode';
    extendedMod.useHybrid = false;
    extendedMod = extendedMod.formPropensitiesGeneral('ODEEricCyt');

    %%    STEP 3.C. -- Change parameters manually, solve, and make plots
    % extendedMod.parameters(4,:) = {'knuc2cyt',0.027};
    % extendedMod.parameters(15,:) = {'knucdeg',0.001};
    % extendedMod.parameters(16,:) = {'degCyt',0.003};
    if makePlots
        soln = extendedMod.solve;
        plotODEresults(extendedMod,soln,ModelGRfit{dexModel})
    end
    %%    STEP 3.D.1. -- Fit new parameters to match all ODE data.
    % extendedMod.fittingOptions.modelVarsToFit = [1:4,14:18];
    extendedMod.fittingOptions.modelVarsToFit = [4,14:15];
    extendedMod.fittingOptions.logPrior = ...
        @(x)-sum((log10(x(:))'-log10PriorMean(extendedMod.fittingOptions.modelVarsToFit)).^2./...
        (2*log10PriorStd(extendedMod.fittingOptions.modelVarsToFit).^2));

    pars = [extendedMod.parameters{extendedMod.fittingOptions.modelVarsToFit,2}];
    % load('BestParsODEonly062725','pars');
    for i = 1:fitIters
        pars = extendedMod.maximizeLikelihood(pars,fitOptions);
    end
    % Best fit so far: 386.17, but distributions and parameters are terrible.

    %    STEP 3.D.2. -- Plot ODE fit results.
    extendedMod.parameters(extendedMod.fittingOptions.modelVarsToFit,2) = num2cell(pars);
    if makePlots
        soln = extendedMod.solve;
        plotODEresults(extendedMod,soln,ModelGRfit{3})
    end
    fullPars = [extendedMod.parameters{[1:12,14:16],2}];
    %%    STEP 3.F.1. -- Create New Objective Function Combining all of the previous ones.
    % Remove all priors from individual models.
    ModelGRDusp100nM_ext_red.fittingOptions.logPrior = [];
    extendedMod.fittingOptions.logPrior = [];
    % fullPars = [extendedMod.parameters{[1:12,14:18],2}];
    % load('fullPars_062725','fullPars')
    % otherewise, we will use the set that was saved in the data dump from
    % the older workspace.

end
%%    STEP 3.F.2.A -- Fit all objective functions at once -- RNA ONLY
% First, start with just the early time DUSP1 and the ODE model.
logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:4,14,15])).^2./(2*log10PriorStd([1:4,14,15]).^2));

% extendedMod.fittingOptions.modelVarsToFit = [1:12,14];
% Organization = {ModelGRfit{1},[5:12],[5:12],'computeLikelihood',1;...
%     ModelGRfit{2},[5:12],[5:12],'computeLikelihood',1;...
%     ModelGRfit{3},[5:12],[5:12],'computeLikelihood',1;...
%     ModelGRDusp100nM_ext_red,[1:12,14],[1:13],'computeLikelihood',1;...
%     extendedMod,[1:12,14:18],[1:17],'computeLikelihoodODE',0.01};
% load('fullPars_062725','fullPars')
fullPars(5:12)=GRpars; % Update to make sure consistent with GR fit.
extendedMod.parameters(1:12,2) = num2cell(fullPars(1:12));
extendedMod.parameters(14:15,2) = num2cell(fullPars(13:14));
extendedMod.parameters(16:20,2) = {1e-9;1;1;1;1e-9};

Organization = {ModelGRDusp100nM_ext_red,[1:4,14],[1:5],'computeLikelihood',1;...
    extendedMod,[1:4,14,15],[1:6],'computeLikelihoodODE',1};

% noProtein= fullPars([1:4,13:14]).*(1+0.1*randn(size(fullPars([1:4,13:14]))));
noProtein= fullPars([1:4,13:14]);
% Run getTotalFitErr to initialize:
Organization = getTotalFitErr(Organization,noProtein,true);
% Run getTotalFitErr again to calculate initial error:
getTotalFitErr(Organization,noProtein,false)

objAll = @(x)-getTotalFitErr(Organization,exp(x),false)-logPriorAll(exp(x));

for i=1:fitIters
    noProteinLog = log(max(1e-9,noProtein)); noProtein = exp(fminsearch(objAll,noProteinLog,fitOptions));
end
fullPars([1:4,13:17]) = [noProtein,1e-9,1,1e-9];
% Best error so far: 4238.92 (with GRpars fixed from above)
% noProtein=
% 0.428604195715366
% 0.00970870640233206
% 6.80295712014968
% 0.065344009835494
% 2.78527860976539e-05
% 0.0308589791156274
%%    STEP 3.F.2.A -- Fit all objective functions at once -- With Protein Feedback
% First, start with just the early time DUSP1 and the ODE model.

% extendedMod.fittingOptions.modelVarsToFit = [1:12,14];
% Organization = {ModelGRfit{1},[5:12],[5:12],'computeLikelihood',1;...
%     ModelGRfit{2},[5:12],[5:12],'computeLikelihood',1;...
%     ModelGRfit{3},[5:12],[5:12],'computeLikelihood',1;...
%     ModelGRDusp100nM_ext_red,[1:12,14],[1:13],'computeLikelihood',1;...
%     extendedMod,[1:12,14:18],[1:17],'computeLikelihoodODE',0.01};

% logPriorAll = @(x)-sum((log10(x)-log10PriorMean([4,14:18])).^2./(2*log10PriorStd([4,14:18]).^2));
% Organization = {ModelGRDusp100nM_ext_red,[4,14],[1:2],'computeLikelihood',1;...
%     extendedMod,[4,14:18],[1:6],'computeLikelihoodODE',1};
% withProtein = fullPars([4,13:17]); 
% 
% logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:4,14:18])).^2./(2*log10PriorStd([1:4,14:18]).^2));
% Organization = {ModelGRDusp100nM_ext_red,[1:4,14],[1:5],'computeLikelihood',1;...
%     extendedMod,[1:4,14:18],[1:9],'computeLikelihoodODE',1};
% withProtein = fullPars([1:4,13:17]); 

% logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:4,14:20])).^2./(2*log10PriorStd([1:4,14:20]).^2));
% Organization = {ModelGRDusp100nM_ext_red,[1:4,14],[1:5],'computeLikelihood',1;...
%     extendedMod,[1:4,14:20],[1:11],'computeLikelihoodODE',1};
% withProtein = fullPars([1:4,13:19]); 
% 
logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:4,14:16])).^2./(2*log10PriorStd([1:4,14:16]).^2));
Organization = {ModelGRDusp100nM_ext_red,[1:4,14],[1:5],'computeLikelihood',1;...
    extendedMod,[1:4,14:16],[1:7],'computeLikelihoodODE',1};
withProtein = fullPars([1:4,13:15]); 
% 
Organization = getTotalFitErr(Organization,withProtein,true);
getTotalFitErr(Organization,withProtein,false)

objAll = @(x)-getTotalFitErr(Organization,exp(x),false)-logPriorAll(exp(x));

for i=1:fitIters
    withProteinLog = log(max(1e-9,withProtein)); withProtein = exp(fminsearch(objAll,withProteinLog,fitOptions));
end
% Best error so far: 112525 (with GRpars fixed from above)
% fullPars([4,13:17]) = withProtein;
% fullPars([1:4,13:17]) = withProtein;
% fullPars([1:4,13:19]) = withProtein;

%%
% logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:12,14:18])).^2./(2*log10PriorStd([1:12,14:19]).^2));
% 
% Organization = {combinedGRModel.SSITModels{1},[5:12],[5:12],'computeLikelihood',1;...
%     combinedGRModel.SSITModels{2},[5:12],[5:12],'computeLikelihood',1;...
%     combinedGRModel.SSITModels{3},[5:12],[5:12],'computeLikelihood',1;...
%     ModelGRDusp100nM_ext_red,[1:12,14],[1:13],'computeLikelihood',1;...
%     extendedMod,[1:12,14:18],[1:17],'computeLikelihoodODE',1};
% 
% Organization = getTotalFitErr(Organization,fullPars,true);
% getTotalFitErr(Organization,fullPars,false)
% 
% objAll = @(x)-getTotalFitErr(Organization,exp(x),false)-logPriorAll(exp(x));
% %%
% load('fullPars_062525B');
% fitOptions = optimset('Display','iter','MaxIter',500);
% for i = 1:5
%     fullParsLog = log(max(1e-6,fullPars));
%     fullPars = exp(fminsearch(objAll,fullParsLog,fitOptions));
%     save('fullPars_062525B','fullPars')
% end

%%    STEP 3.F.3. -- Plot all results
% Plot GR Distribution
if makePlots
    figNums = [1:12];
    for i = figNums
        try
            close(i);
        catch
        end
    end
    makeGRPlots(combinedGRModel,fullPars(5:12),GR_Data);
end
%%
% Plot DUSP1 100nm FIT and other PREDICTED distributions
if makePlots
    showCases = [1,1,1,0];
    ModelGRDusp100nM_ext_red.fittingOptions.modelVarsToFit = [1:12,14];
    ModelGRDusp100nM_ext_red.parameters([1:12,14],2) = num2cell(fullPars([1:13]));
    makePlotsDUSP1({ModelGRDusp100nM_ext_red},ModelGRDusp100nM_ext_red,fullPars([1:13]),Dusp1FitCases,showCases,Dusp1Data)
end
%%
% Plot ODE Cyt FIT Results at 100nM Dex
% extendedMod.parameters([1:12,14:18],2) = num2cell(fullPars);
extendedMod.parameters([1:12,14:16],2) = num2cell(fullPars);

soln100 = extendedMod.solve;
if makePlots
    plotODEresults(extendedMod,soln100,ModelGRfit{3},500)
    set(gcf,'Name','ODE Fits -- 100nM Dex')
end

% Plot ODE Predictions at other DEX concentrations
extendedMod0p3 = extendedMod.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],['(TAB.dex_conc==0.3|TAB.dex_conc==0)']});
extendedMod0p3.parameters(13,:) = {'Dex0',0.3};

extendedMod1 = extendedMod.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],['(TAB.dex_conc==1|TAB.dex_conc==0)']});
extendedMod1.parameters(13,:) = {'Dex0',1.0};

extendedMod10 = extendedMod.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],['(TAB.dex_conc==10|TAB.dex_conc==0)']});
extendedMod10.parameters(13,:) = {'Dex0',10};

if makePlots
    plotODEresults(extendedMod1,extendedMod1.solve,ModelGRfit{1},501)
    set(gcf,'Name','ODE Predictions -- 1.0nM Dex')

    plotODEresults(extendedMod10,extendedMod10.solve,ModelGRfit{2},502)
    set(gcf,'Name','ODE Predictions -- 10nM Dex')
end
%%    STEP 3.G.1. -- Create SSA model to predict cytoplasmic distributions
SSAModel_100 = extendedMod;
SSAModel_100.solutionScheme = 'SSA';
SSAModel_100 = SSAModel_100.formPropensitiesGeneral('EricSSA');
SSAModel_100.ssaOptions.verbose = false;
% SSAModel_100.tSpan = [-500,ModelGRDusp100nM.tSpan];
SSAModel_100.tSpan = [-500,0:10:500];
% A negative initial time is needed to allow model to equilibrate before
% starting.  This causes long run times.
% Note that because the SSA creates independent sets for each timepoint, I
% am doin this to constrain it to just 100 simulations total.
% SSAModel_100.initialCondition = [2;0;round(soln100.ode(1,3:7))'];
SSAModel_100.initialCondition = [2;0;round(soln100.ode(1,3:6))'];
SSAModel_100.initialTime = SSAModel_100.tSpan(1);
SSAModel_100.hybridOptions = [];
SSAModel_100.useHybrid = false;
SSAModel_100.ssaOptions.useParallel = true;
SSAModel_100 = SSAModel_100.formPropensitiesGeneral('SSAExtendedModel');
%    STEP 3.G.2. -- Run SSA Simulations
SSAModel_100.tSpan = [-500,extendedMod.tSpan];
SSAModel_100.ssaOptions.nSimsPerExpt = 1000/length(SSAModel_100.tSpan);

%%
indsPars = [1:4,15:16];
x =[SSAModel_100.parameters{indsPars,2}];
objCyt = @(x)computeCytError(exp(x),indsPars,SSAModel_100,extendedMod);
objCyt(log(x))

for i=1:3
    x = exp(fminsearch(objCyt,log(x),fitOptions));
end

SSAModel_100.parameters(indsPars,2) = num2cell(x);

%%
ssaSoln_100 =SSAModel_100.solve;
%    STEP 3.G.3. -- Plot SSA Results for Cytoplasmic Distributions (100nM Dex)
% if makePlots
    % Plot Fits for 100nM Dex
    f=figure(600); clf; set(f,'Name','Nuclear RNA, 100nM Dex')
    makeCytDistPlots(ssaSoln_100,extendedMod,600,[1:12],5,1,[0:5:300],true)
    f=figure(601); clf; set(f,'Name','Cytoplasmic RNA, 100nM Dex')
    sum(makeCytDistPlots(ssaSoln_100,extendedMod,601,[1:12],6,2,[0:5:300],true))
    f=figure(602); clf; set(f,'Name','Nuc vs. Cyto RNA, 100nM Dex')
    makeNucCytScatterPlots(ssaSoln_100,extendedMod,602,[1:12],[5,6],[1,2],true)
% end

%%    STEP 3.F.1. -- Predict Cyt distributions for other 0.3nM Dex 
SSAModel_0p3 = SSAModel_100;
SSAModel_0p3.parameters(13,:) = {'Dex0',0.3};
ssaSoln_0p3 = SSAModel_0p3.solve;

%    STEP 3.F.2. -- Predict Cyt distributions for other 1.0nM Dex 
SSAModel_1 = SSAModel_100;
SSAModel_1.parameters(13,:) = {'Dex0',1.0};
ssaSoln_1 = SSAModel_1.solve;

%    STEP 3.F.3. -- Predict Cyt distributions for other 10nM Dex 
SSAModel_10 = SSAModel_100;
SSAModel_10.parameters(13,:) = {'Dex0',10};
ssaSoln_10 = SSAModel_10.solve;
% if makePlots
    %    STEP 3.F.4. -- Make resulting plots
    f=figure(603); clf; set(f,'Name','Nuclear RNA, 0.3nM Dex')
    f=figure(605); clf; set(f,'Name','Cytoplasmic RNA, 0.3nM Dex')
    f=figure(605); clf; set(f,'Name','Nuc vs. Cyto RNA, 0.3nM Dex')
    makeCytDistPlots(ssaSoln_0p3,extendedMod0p3,603,[1:7],5,1,[0:15:300],true)
    makeCytDistPlots(ssaSoln_0p3,extendedMod0p3,604,[1:7],6,2,[0:15:300],true)
    makeNucCytScatterPlots(ssaSoln_0p3,extendedMod0p3,605,[1:7],[5,6],[1,2],true)

    f=figure(606); clf; set(f,'Name','Nuclear RNA, 1.0nM Dex')
    f=figure(607); clf; set(f,'Name','Cytoplasmic RNA, 1.0nM Dex')
    f=figure(608); clf; set(f,'Name','Nuc vs. Cyto RNA, 1.0nM Dex')
    makeCytDistPlots(ssaSoln_1,extendedMod1,606,[1:7],5,1,[0:15:300],true)
    makeCytDistPlots(ssaSoln_1,extendedMod1,607,[1:7],6,2,[0:15:300],true)
    makeNucCytScatterPlots(ssaSoln_1,extendedMod1,608,[1:7],[5,6],[1,2],true)

    f=figure(609); clf; set(f,'Name','Nuclear RNA, 10nM Dex')
    f=figure(610); clf; set(f,'Name','Cytoplasmic RNA, 10nM Dex')
    f=figure(611); clf; set(f,'Name','Nuc vs. Cyto RNA, 10nM Dex')
    makeCytDistPlots(ssaSoln_10,extendedMod10,609,[1:7],5,1,[0:15:300],true)
    makeCytDistPlots(ssaSoln_10,extendedMod10,610,[1:7],6,2,[0:15:300],true)
    makeNucCytScatterPlots(ssaSoln_10,extendedMod10,611,[1:7],[5,6],[1,2],true)
% end

% %%
% %% Save Results for Easier Use in subsequent runs.
% %parsAll_GR_Dusp1_TS = [extendedMod.parameters{:,2}];
% %parsAll_GR_Dusp1_TS(16) = parsAllandTS(end);
% varNames = unique({'ModelGR'
%     'GRfitCases'
%     'log10PriorMean'
%     'log10PriorStd'
%     'GRpars'
%     'ModelGRparameterMap'
%     'ModelGRfit'
%     'boundGuesses'
%     'ModelGRDusp100nM'
%     'GRfitCases'
%     'log10PriorMean'
%     'log10PriorStd'
%     'duspLogPrior'
%     'DUSP1pars'
%     'Dusp1FitCases'
%     'ModelGRfit'
%     'extendedMod'
%     'ModelGRDusp100nM_ext_red'
%     'fullPars'
%     'fimResults'
%     %'DUSP1parsIntensity'
%     });
% 
% save('workspaceJuly24',varNames{:})
% 
% %% Extra Function
% 
