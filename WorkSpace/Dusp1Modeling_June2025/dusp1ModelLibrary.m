function [Model,log10PriorMean,log10PriorStd] = dusp1ModelLibrary(modelName,regenerate,modelLibrary)
arguments
    modelName = 'ModelBase'
    regenerate = false
    modelLibrary = 'savedParameters/GRDusp1ModelLibrary'
end

if regenerate

    GR_Data = 'RonData062025/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv';
    Dusp1Data = 'RonData062025/ssit_gated_Jun19';    
    
    ModelBase = SSIT;
    ModelBase.species = {
        'offGene';
        'onGene' ;
        'cytGR'  ;
        'nucGR'  ;
        'rna'    ;
        'rCyt'   };
    ModelBase.initialCondition = [2;0;24;1;5;0];

    ModelBase.propensityFunctions = {
        '(kon*nucGR)*offGene';
        'koff*onGene';
        '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';
        'knc*nucGR';
        'kg1';
        'gg1*cytGR';
        'gg2*nucGR';
        'kr*onGene';
        'knuc2cyt*rna';
        '(degCyt0 + degCyt1./(1+degCytA*rCyt) ).*rCyt'};
    ModelBase.stoichiometry = [
     -1     1     0     0     0     0     0     0     0    0;
     1    -1     0     0     0     0     0     0     0     0;
     0     0    -1     1     1    -1     0     0     0     0;
     0     0     1    -1     0     0    -1     0     0     0;
     0     0     0     0     0     0     0     1    -1     0;
     0     0     0     0     0     0     0     0     1    -1];

    ModelBase.parameters = {'koff',0.928272906611829;
        'kon',0.00435679746751789;
        'kr',13.6899693803118;
        'knuc2cyt',0.0303917515033354;
        'kcn0',0.00605936338887321;
        'kcn1',0.100801493996861;
        'gDex',3.06931711147851e-05;
        'knc',0.0101878879366104;
        'kg1',0.109965006106303;
        'gg1',0.0035378622695277;
        'gg2',0.00955287461619338;
        'MDex',14.1807103753798;
        'Dex0',100;
        'degCyt0',0.005;
        'degCyt1',0.005;
        'degCytA',0.01};

    % Priors (based on magnitude of previously guessed values).
    log10PriorMean = [
        0, ...'koff',0.928272906611829;
        -3, ...'kon',0.00435679746751789;
        1, ...'kr',13.6899693803118;
        -2, ...'knuc2cyt',0.0303917515033354;
        -2, ...'kcn0',0.00605936338887321;
        -1, ...'kcn1',0.100801493996861;
        -5, ...'gDex',3.06931711147851e-05;
        -2, ...'knc',0.0101878879366104;
        -1, ...'kg1',0.109965006106303;
        -3, ...'gg1',0.0035378622695277;
        -2, ...'gg2',0.00955287461619338;
        1, ...'MDex',14.1807103753798;
        NaN, ...'Dex0',100;
        -2, ...'degCyt0',0.00834750349323117;
        -2, ... %'degCyt1',0.00927999553405981
        -2]';    %'degCytA'
    log10PriorStd = 2*ones(16,1);

    ModelBase.fspOptions.initApproxSS = true;
    ModelBase.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
    ModelBase.sensOptions.solutionMethod = 'finiteDifference';
 
    %% Simplify to GR model
    ModelGR = ModelBase;
    GRspeciesInds = [3,4];
    GRreactionInds = (3:7);
    GRparameterInds = (5:12);
    ModelGR.species = ModelGR.species(GRspeciesInds);
    ModelGR.initialCondition = ModelGR.initialCondition(GRspeciesInds);
    ModelGR.propensityFunctions = ModelGR.propensityFunctions(GRreactionInds);
    ModelGR.stoichiometry = ModelGR.stoichiometry(GRspeciesInds,GRreactionInds);
    ModelGR.fittingOptions.modelVarsToFit = (5:12);
    ModelGR.customConstraintFuns = {'cytGR+nucGR'};
    ModelGR = ModelGR.formPropensitiesGeneral('ModelGR');
    ModelGR.fspOptions.bounds = [0;0;30;30;30];
    [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve;
    [~,ModelGR.fspOptions.bounds] = ModelGR.solve(FSPGrSoln.stateSpace);

    %% Add GR Data and create multimodel
    GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
        '10','10',102,'GR Fit (10nM Dex)';...
        '100','100',103,'GR Fit (100nM Dex)'};
    ModelGRparameterMap = cell(1,size(GRfitCases,1));
    ModelGRfit = cell(1,size(GRfitCases,1));
    boundGuesses = cell(1,size(GRfitCases,1));
    for i=1:3
        if i==3 % Include zero time data within the 100nM dataset. Exclude times (20,40,60,90,150) due to missing replicas.
            ModelGRfit{i} = ModelGR.loadData(GR_Data,...
                {'cytGR','normGRcyt';'nucGR','normGRnuc'},...
                {[],[], ...
                ['(TAB.dex_conc==',GRfitCases{i,1},'|(TAB.dex_conc==0&TAB.time==0))&TAB.time~=20&TAB.time~=40&TAB.time~=60&TAB.time~=90&TAB.time~=150']});
        else % Exclude zero time data within the other datasets to avoid double counting in fitting. 
             % Exclude times (20,40,60,90,150) due to missing replicas.
            ModelGRfit{i} = ModelGR.loadData(GR_Data,...
                {'cytGR','normGRcyt';'nucGR','normGRnuc'},...
                {[],[], ...
                ['TAB.dex_conc==',GRfitCases{i,1},'&TAB.time~=20&TAB.time~=40&TAB.time~=60&TAB.time~=90&TAB.time~=150']});
        end

        ModelGRfit{i}.parameters(13,:) = {'Dex0', str2num(GRfitCases{i,1})};
        ModelGRfit{i}.tSpan = ModelGRfit{i}.dataSet.times;
        ModelGRparameterMap(i) = {(1:8)};
        boundGuesses{i} = [0;0;30;30;30];
    end

    logPriorGR = @(x)-sum((log10(reshape(x,[numel(x),1]))-log10PriorMean(GRparameterInds)).^2./(2*log10PriorStd(GRparameterInds).^2));
    combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap,logPriorGR);
    combinedGRModel = combinedGRModel.initializeStateSpaces(boundGuesses);

    %%  Nuclear DUSP1 Model
    ModelDUSP1 = ModelBase;
    DUSP1speciesInds = (1:5);
    DUSP1reactionInds = (1:9);
    DUSP1parameterInds = (1:12);
    ModelDUSP1.species = ModelDUSP1.species(DUSP1speciesInds);
    ModelDUSP1.initialCondition = ModelDUSP1.initialCondition(DUSP1speciesInds);
    ModelDUSP1.propensityFunctions = ModelDUSP1.propensityFunctions(DUSP1reactionInds);
    ModelDUSP1.stoichiometry = ModelDUSP1.stoichiometry(DUSP1speciesInds,DUSP1reactionInds);
    ModelDUSP1.fittingOptions.modelVarsToFit = (1:4);
    ModelDUSP1.useHybrid = true;
    ModelDUSP1.hybridOptions.upstreamODEs = {'cytGR','nucGR'};
    ModelDUSP1 = ModelDUSP1.formPropensitiesGeneral('ModelDUSP1');
    ModelDUSP1.fspOptions.bounds = [0;0;0;2;2;500];
    
    [FSPDusp1Soln,ModelDUSP1.fspOptions.bounds] = ModelDUSP1.solve;
    [~,ModelDUSP1.fspOptions.bounds] = ModelDUSP1.solve(FSPDusp1Soln.stateSpace);    
    
    ModelDUSP1.fittingOptions.logPrior = ...
        @(x)-sum((log10(reshape(x,[numel(x),1]))- log10PriorMean(ModelDUSP1.fittingOptions.modelVarsToFit)).^2./ ...
        (2*log10PriorStd(ModelDUSP1.fittingOptions.modelVarsToFit).^2));
    
    ModelDUSP1_100nM = ModelDUSP1.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'},...
        {[],[],['(TAB.dex_conc==100|(TAB.dex_conc==0&TAB.time==0))', ...
        '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});
    ModelDUSP1_100nM.parameters{13,2} = 100;

    ModelDUSP1_10nM = ModelDUSP1.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'},...
        {[],[],['(TAB.dex_conc==10|(TAB.dex_conc==0&TAB.time==0))', ...
        '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});
    ModelDUSP1_10nM.parameters{13,2} = 10;
    
    ModelDUSP1_1nM = ModelDUSP1.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'},...
        {[],[],['(TAB.dex_conc==1|(TAB.dex_conc==0&TAB.time==0))', ...
        '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});
    ModelDUSP1_1nM.parameters{13,2} = 1;
    
    ModelDUSP1_0p3nM = ModelDUSP1.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'},...
        {[],[],['(TAB.dex_conc==0.3|(TAB.dex_conc==0&TAB.time==0))', ...
        '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});
    ModelDUSP1_0p3nM.parameters{13,2} = 0.3;

    DexConc = 10.^[-3,-2,-1,0,1,2,3,4];
    DecConcStr = {'0.001','0.01','0.1','1','10','100','1000','10000'};
    ModelPredDexTtr = cell(size(DexConc));
    for i = 1:length(DexConc)
        ModelPredDexTtr{i} = ModelDUSP1.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'},...
        {[],[],['(TAB.dex_conc==',DecConcStr{i},'&TAB.time==75)', ...
        '&TAB.cyto_area>=12593&TAB.cyto_area<=17685']});
        ModelPredDexTtr{i}.parameters{13,2} = DexConc(i);

    end

    %% Full ODE Model
    fullODEModel = ModelBase;
    fullODEModel.solutionScheme = 'ode';
    fullODEModel = fullODEModel.formPropensitiesGeneral('fullODEModel');

    fullODEModel_0p3 = fullODEModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==0.3|(TAB.dex_conc==0&TAB.time==0))'});
    fullODEModel_0p3.parameters(13,:) = {'Dex0',0.3};

    fullODEModel_1 = fullODEModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==1|(TAB.dex_conc==0&TAB.time==0))'});
    fullODEModel_1.parameters(13,:) = {'Dex0',1.0};

    fullODEModel_10 = fullODEModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==10|(TAB.dex_conc==0&TAB.time==0))'});
    fullODEModel_10.parameters(13,:) = {'Dex0',10};

    fullODEModel_100 = fullODEModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==100|(TAB.dex_conc==0&TAB.time==0))'});
    fullODEModel_100.parameters(13,:) = {'Dex0',100};

    %% Full SSA Model 
    soln100 = fullODEModel_100.solve;
    fullSSAModel = ModelBase;
    fullSSAModel.solutionScheme = 'SSA';
    fullSSAModel.tSpan = [-500,fullSSAModel.tSpan];
    fullSSAModel.initialCondition = [2;0;round(soln100.ode(1,3:6))'];
    fullSSAModel = fullSSAModel.formPropensitiesGeneral('fullSSAModel');
    fullSSAModel.ssaOptions.useTimeVar = true;
    fullSSAModel.ssaOptions.signalUpdateRate = 1;
    fullSSAModel.initialTime = -500;
    fullSSAModel.ssaOptions.useParallel = true;

    fullSSAModel_0p3 = fullSSAModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==0.3|(TAB.dex_conc==0&TAB.time==0))'});
    fullSSAModel_0p3.parameters(13,:) = {'Dex0',0.3};

    fullSSAModel_1 = fullSSAModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==1|(TAB.dex_conc==0&TAB.time==0))'});
    fullSSAModel_1.parameters(13,:) = {'Dex0',1};

    fullSSAModel_10 = fullSSAModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==10|(TAB.dex_conc==0&TAB.time==0))'});
    fullSSAModel_10.parameters(13,:) = {'Dex0',10};
    
    fullSSAModel_100 = fullSSAModel.loadData(Dusp1Data,...
        {'rna','num_nuc_spots'; ...
        'rCyt','num_cyto_spots'},...
        {[],[],'(TAB.dex_conc==100|(TAB.dex_conc==0&TAB.time==0))'});
    fullSSAModel_100.parameters(13,:) = {'Dex0',100};

    %%

    save(modelLibrary,'ModelBase','ModelGR','combinedGRModel','ModelDUSP1_0p3nM',...
        'ModelDUSP1_1nM','ModelDUSP1_10nM','ModelDUSP1_100nM','ModelPredDexTtr',...
        'fullODEModel_0p3','fullODEModel_1','fullODEModel_10','fullODEModel_100',...
        'fullSSAModel_0p3','fullSSAModel_1','fullSSAModel_10','fullSSAModel_100',...
        'log10PriorMean','log10PriorStd')
end
Models = load(modelLibrary,modelName,'log10PriorMean','log10PriorStd');
Model = Models.(modelName);
log10PriorMean = Models.log10PriorMean;
log10PriorStd = Models.log10PriorStd;
