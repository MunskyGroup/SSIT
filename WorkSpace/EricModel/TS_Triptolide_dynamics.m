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
    try
        ModelGRDusp100nM_ext_red.propensitiesGeneral{1}.stateDependentFactor(0);
    catch
        ModelGRDusp100nM_ext_red = ModelGRDusp100nM_ext_red.formPropensitiesGeneral('ExtModel100nm');
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
fitOptions.MaxIter = 10;
fitOptions.UseParallel = true;
for i = 1:1
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
showCases = [1,1,1,1];
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

%%
%%  STEP 5 -- Tryptolide Peturbation
ktptl = log(2)/4;
%%    STEP 5.A. -- DUSP1 Model
ModelGRDusp100nM_ext_red_TPL = ModelGRDusp100nM_ext_red;
ModelGRDusp100nM_ext_red_TPL.propensityFunctions{8} = 'kr*onGene*Itrpt';
% [0,20,75,150,180] % choose the time of the TPT treatment.
tptlTime = '75';
ModelGRDusp100nM_ext_red_TPL.inputExpressions(2,:) = {'Itrpt',['(t<=',tptlTime,') + (t>',tptlTime,')*exp(-ktptl*(t-',tptlTime,'))']};
ModelGRDusp100nM_ext_red_TPL.parameters(15,:) = {'degCyt',NaN};
ModelGRDusp100nM_ext_red_TPL.parameters(16,:) = {'ktptl',ktptl};
ModelGRDusp100nM_ext_red_TPL.dataSet = [];
ModelGRDusp100nM_ext_red_TPL = ModelGRDusp100nM_ext_red_TPL.loadData('EricData/TryptolideData.csv',...
            {'rna','RNA_DUSP1_nuc'},...
            {'Time_TPL',tptlTime});
ModelGRDusp100nM_ext_red_TPL = ModelGRDusp100nM_ext_red_TPL.formPropensitiesGeneral('TrptModel');
ModelGRDusp100nM_ext_red_TPL.makeFitPlot

%%    STEP 5.B. -- Nuc/Cyt ODE Model 
extendedMod_tptl = extendedMod;
extendedMod_tptl.propensityFunctions{8} = 'kr*onGene*Itrpt';
extendedMod_tptl.inputExpressions(2,:) = {'Itrpt',['(t<=',tptlTime,') + (t>',tptlTime,')*exp(-ktptl*(t-',tptlTime,'))']};
extendedMod_tptl.parameters(16,:) = {'ktptl',ktptl};
extendedMod_tptl.dataSet = [];
extendedMod_tptl = extendedMod_tptl.loadData('EricData/TryptolideData.csv',...
            {'rna','RNA_DUSP1_nuc';...
            'rCyt','RNA_DUSP1_cyto'},...
            {'Time_TPL',tptlTime});
extendedMod_tptl = extendedMod_tptl.formPropensitiesGeneral('TrptModelODE');
plotODEresults(extendedMod_tptl,extendedMod_tptl.solve,ModelGRfit{3},501)
%%    STEP 5.C. -- TS after Tryptolide
TSMod_tptl = ModelTS;
TSMod_tptl.propensityFunctions{8} = 'kr*onGene*Itrpt';
TSMod_tptl.inputExpressions(2,:) = {'Itrpt',['(t<=',tptlTime,') + (t>',tptlTime,')*exp(-ktptl*(t-',tptlTime,'))']};
TSMod_tptl.parameters(15,:) = {'degCyt',NaN};
TSMod_tptl.parameters(16,:) = {'ktptl',ktptl};  % Rate of diffusion of TPT to nucleus.
TSMod_tptl = TSMod_tptl.formPropensitiesGeneral('TrptModelODE');

dex = 100;
allData = readtable('EricData/TryptolideData.csv');
allData = allData(allData.Time_TPL==str2double(tptlTime),:);

TSMod_tptl.tSpan = unique(allData.Time_index)';

[TSLikelihood,modelResults,dataResults] = computeTSlikelihood(parsAllandTS(indsTSpars),TSMod_tptl,allData,dex,52,true); 
makeTSPlots(modelResults,dataResults,TSMod_tptl,TSthresh,[1200,1201])

%%    STEP 5.D. -- Add ODE model after TPL to Objective Function
% Remove redundant times from TPL model.
% Currently this doesn not seem to be working.   Need to ignore some time
% poitns.  the folloiwing should give two different results with the
% seconde being smaller.
ModelGRDusp100nM_ext_red_TPL_fit = ModelGRDusp100nM_ext_red_TPL;

% Remove redundant times from the data set (all times before application
% of TPL);
tptlTime = '75';
ModelGRDusp100nM_ext_red_TPL_fit.inputExpressions(2,:) = {'Itrpt',['(t<=',tptlTime,') + (t>',tptlTime,')*exp(-ktptl*(t-',tptlTime,'))']};
[~,J] = intersect(ModelGRDusp100nM_ext_red_TPL_fit.dataSet.times,[85;90;105;135]);
ModelGRDusp100nM_ext_red_TPL_fit.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor = ...
ModelGRDusp100nM_ext_red_TPL.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(J,:);
ModelGRDusp100nM_ext_red_TPL_fit.dataSet.times = ModelGRDusp100nM_ext_red_TPL.dataSet.times(J);
ModelGRDusp100nM_ext_red_TPL_fit.tSpan = unique([0,ModelGRDusp100nM_ext_red_TPL_fit.dataSet.times]);
ModelGRDusp100nM_ext_red_TPL_fit.dataSet.nCells = ModelGRDusp100nM_ext_red_TPL.dataSet.nCells(J);
ModelGRDusp100nM_ext_red_TPL_fit.dataSet.mean = ModelGRDusp100nM_ext_red_TPL.dataSet.mean(J);

OrganizationTPL = {ModelGRfit{1},[5:12],[5:12],'computeLikelihood',1;...
    ModelGRfit{2},[5:12],[5:12],'computeLikelihood',1;...
    ModelGRfit{3},[5:12],[5:12],'computeLikelihood',1;...
    ModelGRDusp100nM_ext_red,indsDuspMod,indsDuspDat,'computeLikelihood',1;...
    ModelGRDusp100nM_ext_red_TPL_fit,[indsDuspMod,16],[indsDuspDat,15],'computeLikelihood',1;...
    extendedMod,indsODEmod,indsODEdat,'computeLikelihoodODE',0.01;...
    };

log10PriorMean(17) = -1; % tpl effectivity rate
log10PriorStd(17) = 2;
logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:12,14:17])).^2./(2*log10PriorStd([1:12,14:17]).^2));

objTPL = @(x)-getTotalFitErr(OrganizationTPL,exp(x),false,extraObjs)-logPriorAll(exp(x));

if loadPrevious
     varNamesTPL = {'parsTPL'};
     load(savedWorkspace,varNamesTPL{:})
else
    parsTPL=parsAllandTS;
    parsTPL(16) = log(2)/5;
end

% Check that the function works.
tic
    objTPL(log(parsTPL))
toc
%%    STEP 5.E. -- Run the fit with fminsearch
fitOptions.MaxIter = 1000;
fitOptions.UseParallel = true;
for i = 1:5
    parsTPL = exp(fminsearch(objTPL,log(parsTPL),fitOptions));
    bestObj = objTPL(log(parsTPL));
end

%%
ModelGRDusp100nM_ext_red_TPL.parameters([indsDuspMod,16],2) = num2cell(parsTPL([indsDuspDat,15]));
% ModelGRDusp100nM_ext_red_TPL.makeFitPlot

showCases = [0,0,0,1];
makePlotsDUSP1({ModelGRDusp100nM_ext_red},ModelGRDusp100nM_ext_red,parsTPL([indsDuspDat]),Dusp1FitCases,showCases,parsTPL([16]))

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
    
    subplot(ceil(NT/3),3,iT)
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
set(gca,'fontsize',16,'ylim',[0,0.5])
ylabel('Fraction with TS')
xlabel('time (min)')
%
%    STEP 4.C.2. -- Add data to plots
figure(figNum(1))
for iT = 1:NT
    subplot(ceil(NT/3),3,iT); hold on
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

%%
%% Save Results for Easier Use in subsequent runs.
parsAll_GR_Dusp1_TS = [extendedMod.parameters{:,2}];
parsAll_GR_Dusp1_TS(16) = parsAllandTS(end);
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
    'parsAll_GR_Dusp1_TS'
    'parsTPL'
    });

save('workspaceJuly24',varNames{:})