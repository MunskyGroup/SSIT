% This script fits/predicts all data from Eric Ron's experimental
% measurement of GR translocation and DUSP1 activation/maturation versus
% time after different levels of dexamethasone. 

close all 
addpath(genpath('../../src'));
addpath('tmpPropensityFunctions/')
modelLibrary = 'savedParameters/GRDusp1ModelLibrary';

fitOptions = optimset('Display','iter','MaxIter',10);
fitIters = 1;
makePlots = true;
computeFIM = false;
runFits = false;
runGRMH = false;
runDusp1MH = false;
runAllMH = false;

%% STEP 1. GR Model
% We start by fitting the GR translocation data.
%%   STEP 1.A -- Load GR MultiModel
[combined_GRModel,log10PriorMean,log10PriorStd] = dusp1ModelLibrary('combinedGRModel',false,modelLibrary);
load('savedParameters/GRparameters','GRpars')
%%   STEP 1.B -- Fit GR MultiModel
if runFits
    for jj = 1:fitIters
        combined_GRModel = combined_GRModel.updateModels(GRpars, false);
        GRpars = combined_GRModel.maximizeLikelihood(GRpars, fitOptions);
        save('EricModel_MMDex','GRpars')
    end
end
combined_GRModel = combined_GRModel.updateModels(GRpars, false);
% Best error so far: 33776.3
%%   STEP 1.C -- Plot GR Fit Results
if makePlots
    figNums = [1:12];
    for i = figNums
        try
            close(i);
        catch
        end
    end
    makeGRPlots(combined_GRModel,GRpars,false)
end
%%   STEP 1.D -- Compute FIM for GR parameters (optional)
if computeFIM
    log10PriorStd = 2*ones(1,15);
    combined_GRModel = combined_GRModel.computeFIMs([],'log');
    GRParsInds = combined_GRModel.SSITModels{1}.fittingOptions.modelVarsToFit;
    fimGR_withPrior = combined_GRModel.FIM.totalFIM+... % the FIM in log space.
        diag(1./(log10PriorStd(GRParsInds)*log(10)).^2);  % Add prior in log space.
end

%%   STEP 1.E -- Run MH on GR Models (optional)
% TODO -- not sure that this is working or if the model is just
% underdetermined.
if runGRMH
    MHFitOptions.thin=1;
    MHFitOptions.numberOfSamples=100;
    MHFitOptions.burnIn=0;
    MHFitOptions.progress=true;
    MHFitOptions.numChains = 1;
    MHFitOptions.useFIMforMetHast = true;
    MHFitOptions.saveFile = 'TMPEricMHGR.mat';
    [~,~,MHResultsGR] = combinedGRModel.maximizeLikelihood(...
        GRpars, MHFitOptions, 'MetropolisHastings');
    save('savedParameters/MHResultsGR','MHResultsGR')
    delete(MHFitOptions.saveFile)
end
 
%%
%% STEP 2 -- Nuclear DUSP1 Model
% Next, we consider the nuclear DUSP1 data.
%%   STEP 2.A -- Load Nuclear DUSP1 Model
ModelDUSP1_100nM = dusp1ModelLibrary('ModelDUSP1_100nM',false,modelLibrary);
load('savedParameters/DUSP1parameters','DUSP1pars')
ModelDUSP1_100nM.parameters(ModelDUSP1_100nM.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
%%   STEP 2.B -- Fit Nuclear DUSP1 Model
if runFits
    for i = 1:fitIters
        fitOptions.suppressFSPExpansion = false;
        DUSP1pars = ModelDUSP1_100nM.maximizeLikelihood(DUSP1pars, fitOptions);
        ModelDUSP1_100nM.parameters(ModelDUSP1_100nM.fittingOptions.modelVarsToFit,2) = ...
            num2cell(DUSP1pars);
        save('savedParameters/DUSP1parameters','DUSP1pars')
    end
end
% Best error so far: 17973.9

%%   STEP 2.C -- Plot Nuclear DUSP1 Fits and Predictions
if makePlots
    figNums = [201,301,302,303,221,321,322,323];
    for i = figNums
        try
            close(i);
        catch
        end
    end
    showCases = [1,1,1,0];
    makePlotsDUSP1Simplified(ModelDUSP1_100nM,DUSP1pars,showCases,modelLibrary);
end
%%   STEP 2.D -- Compute FIM for Dusp1 Parameters
if computeFIM||runDusp1MH
    % STEP 2.D.1. -- Compute sensitivity of the FSP solution
    ModelDUSP1_100nM.solutionScheme = 'fspSens';
    sensSoln = ModelDUSP1_100nM.solve();
    ModelDUSP1_100nM.solutionScheme = 'FSP';
    % STEP 2.D.2. -- Compute FIM
    fimResults = ModelDUSP1_100nM.computeFIM(sensSoln.sens,'log');

    % In the following, the log-prior is used as a prior co-variance matrix.
    % This will be used in the FIM calculation as an FIM without new evidence
    % being set equal to the inverse of this covariance matrix.  More rigorous
    % justification is needed to support this heuristic.
    log10PriorStd = 2*ones(1,16);
    fimTotal = ModelDUSP1_100nM.evaluateExperiment(fimResults,ModelDUSP1_100nM.dataSet.nCells,...
        diag(log10PriorStd.^2));

    % Select the portion of the FIM relating to this stage of the model.
    FIMfree = fimTotal{1}(1:4,1:4);
    if min(eig(FIMfree))<1
        disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
        FIMfree = FIMfree + 1*eye(length(FIMfree));
    end
    covFree = FIMfree^-1;
    covFree = 0.5*(covFree+covFree');
end
%
%%   STEP 2.E -- Run Metropolis Hastings Search
% % Run MH.  Will need to increase the number of samples for final
% % figure generation.
if runDusp1MH
    MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
    MHFitOptions.thin=1;
    MHFitOptions.numberOfSamples=1000;
    MHFitOptions.burnIn=0;
    MHFitOptions.progress=true;
    MHFitOptions.numChains = 1;
    MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
    [DUSP1pars,~,MHResultsDusp1] = ModelDUSP1_100nM.maximizeLikelihood(...
        [], MHFitOptions, 'MetropolisHastings');
    save('savedParameters/MHResultsDusp1','MHResultsDusp1')
    delete('TMPEricMHDusp1.mat')
    ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);

    %  STEP 2.F -- Plot the MH results
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
end
%%
%%  STEP 3 -- Nuclear and Cytoplasmic DUSP1 Model
%%   STEP 3.A -- Load SSA Model
% Load and associate parameters into SSA model
SSAModelAll_100 = dusp1ModelLibrary('fullSSAModel_100',false,modelLibrary);
indsAll = [1:12,14:16];
allParsSaveFileName = 'savedParameters/AllParameters';
load(allParsSaveFileName,'bestx');
SSAModelAll_100.parameters(indsAll,2) = num2cell(bestx);
SSAModelAll_100.fittingOptions.modelVarsToFit = indsAll;
% Set number of simulations
SSAModelAll_100.ssaOptions.nSimsPerExpt = 1000/length(SSAModelAll_100.tSpan);

%% Fit Full Model Using SSA and MH
if runAllMH
    % Define objective function
    objAll = @(x)-(computeCytError(x,indsAll,SSAModelAll_100,true,allParsSaveFileName, ...
        true,combined_GRModel.SSITModels{3})+...
        sum((log10(reshape(x,[numel(x),1]))-log10PriorMean(indsAll)).^2./(2*log10PriorStd(indsAll))));

    MHFitOptions.obj = objAll;
    MHFitOptions.logForm = true;
    MHFitOptions.proposalDistribution=@(x)x+0.05*randn(size(x));
    MHFitOptions.saveFile = 'TMPMHll.mat';
    MHFitOptions.numberOfSamples=5000;
    load(allParsSaveFileName,'bestx')
    [~,~,MHResultsAll] = SSAModelAll_100.maximizeLikelihood(bestxAll',MHFitOptions,'MetropolisHastings');
    save('savedParameters/MHResultsAll',MHResultsAll)
    delete(MHFitOptions.saveFile)
end

%% Make Plots using Best Fit Parameter Set
if makePlots
    for iModel = 1:4
        fullModelNames = {'fullSSAModel_0p3','fullSSAModel_1','fullSSAModel_10','fullSSAModel_100'};
        load(allParsSaveFileName,'bestx');
        fullSSAModel = dusp1ModelLibrary(fullModelNames{iModel},false,modelLibrary);
        fullSSAModel.parameters(indsAll,2) = num2cell(bestx);
        fullSSAModel.ssaOptions.nSimsPerExpt = 1000/length(fullSSAModel.tSpan);

        % Solve and make plots.
        ssaSoln =fullSSAModel.solve;
        f=figure(600+10*iModel); clf; set(f,'Name','Nuclear RNA')
        makeCytDistPlots(ssaSoln,fullSSAModel,600+10*iModel,5,1,[0:5:300],true);
        f=figure(601+10*iModel); clf; set(f,'Name','Cytoplasmic RNA')
        makeCytDistPlots(ssaSoln,fullSSAModel,601+10*iModel,6,2,[0:5:300],true);
        f=figure(602+10*iModel); clf; set(f,'Name','Nuc vs. Cyto RNA')
        makeNucCytScatterPlots(ssaSoln,fullSSAModel,602+10*iModel,[5,6],[1,2],true);
    end
end


return
% %%   STEP 3.C -- Fit Cyto Parameters using SSA Model
% ssaModel = 'fullSSAModel_100'; % Choose from 'fullSSAModel_0p3','fullSSAModel_1','fullSSAModel_10','fullSSAModel_100'
% SSAModel_100 = dusp1ModelLibrary(ssaModel,false,'GRDusp1ModelLibrary');
% SSAModel_100.parameters(5:12,2) = num2cell(GRpars);
% SSAModel_100.parameters(1:4,2) = num2cell(DUSP1pars);
% SSAModel_100.parameters(indsCytPars,2) = num2cell(bestx);
% %%
% SSAModel_100.tSpan = unique([-500,fullSSAModel.tSpan]);
% SSAModel_100.ssaOptions.useParallel = true;
% SSAModel_100.ssaOptions.verbose = false;
% 
% SSAModel_100.ssaOptions.nSimsPerExpt = ceil(1000/length(SSAModel_100.tSpan));
% indsCytPars = [1:4,14:16];
% x =[SSAModel_100.parameters{indsCytPars,2}];
% 
% MHFitOptions.thin=1;
% MHFitOptions.numberOfSamples=1000;
% MHFitOptions.burnIn=0;
% MHFitOptions.progress=true;
% MHFitOptions.numChains = 1;
% MHFitOptions.useFIMforMetHast = false;
% MHFitOptions.saveFile = 'TMPEricNucCyt.mat';
% MHFitOptions.proposalDistribution=@(x)x+0.05*randn(size(x));
% 
% objCyt = @(x)-(computeCytError(x,indsCytPars,SSAModel_100,true,'NucCytParsJuly4')+...
%     sum((log10(reshape(x,[numel(x),1]))-log10PriorMean(indsCytPars)).^2./(2*log10PriorStd(indsCytPars))));
% objCyt(bestx)
% 
% MHFitOptions.obj = objCyt;
% MHFitOptions.logForm = true;
%% Run MH to improve fit with SSA.
% SSAModel_100.fittingOptions.modelVarsToFit = indsCytPars;
% load('NucCytParsJuly4','bestx');
% [~,~,MHResultsNucCyt] = SSAModel_100.maximizeLikelihood(bestx',MHFitOptions,'MetropolisHastings');
% delete(MHFitOptions.saveFile);
% load('NucCytParsJuly4','bestx');
% SSAModel_100.parameters(indsCytPars,2) = num2cell(x);
%%     STEP 3.C.1 --  Plot the results
% SSAModel_100.ssaOptions.nSimsPerExpt = 1000/length(SSAModel_100.tSpan);
% [CytError,ssaSoln] = computeCytError(x,indsCytPars,SSAModel_100);
% CytError
% 
% f=figure(600); clf; set(f,'Name','Nuclear RNA')
% makeCytDistPlots(ssaSoln,SSAModel_100,600,5,1,[0:5:300],true);
% f=figure(601); clf; set(f,'Name','Cytoplasmic RNA')
% makeCytDistPlots(ssaSoln,SSAModel_100,601,6,2,[0:5:300],true,'cdf');
% f=figure(602); clf; set(f,'Name','Nuc vs. Cyto RNA')
% makeNucCytScatterPlots(ssaSoln,SSAModel_100,602,[5,6],[1,2],true);


%%


%% EXTRAS






return
%%   STEP 3.D -- FSP Model for Nuc and Cyt
fspModel = dusp1ModelLibrary('fullSSAModel_100',false,'GRDusp1ModelLibrary');
fspModel.tSpan = fspModel.dataSet.times;
fspModel.initialTime = 0;
fspModel.solutionScheme = 'FSP';
fspModel.useHybrid = true;
fspModel.hybridOptions.upstreamODEs = {'cytGR','nucGR'};
fspModel.customConstraintFuns = {'abs(rna-rCyt)'};
fspModel.fspOptions.verbose = true;
fspModel.fspOptions.bounds = [0,0,0,0,2,2,300,300,100]';
fspModel = fspModel.formPropensitiesGeneral('fspNucCyt');
% [~,~,fspModel] = fspModel.solve;


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
    makeGRPlots(combined_GRModel,fullPars(5:12),GR_Data);
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
