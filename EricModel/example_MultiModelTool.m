%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations
close all 
% clear all
addpath('../CommandLine');

%% Create Base Model for GR Only
ModelGR = SSIT;
ModelGR.species = {'cytGR';'nucGR'};
ModelGR.initialCondition = [20;1];

ModelGR.propensityFunctions = {'(kcn0+kcn1*IDex)*cytGR';'knc*nucGR';...
    'kg1';'gg1*cytGR';'gg2*nucGR'};
ModelGR.stoichiometry = [-1,1,1,-1,0;...
                         1,-1,0,0,-1];
ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
    'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',14e-5;...
    'gg1',1e-5;'gg2',1e-6});

ModelGR.fspOptions.initApproxSS = true;

ModelGR.fittingOptions.modelVarsToFit = (5:11);

ModelGR.inputExpressions = {'IDex','100*exp(-gDex*t)*(t>0)'};
[~,ModelGR.fspOptions.bounds] = ModelGR.solve;

%% Load previously fit parameter values (optional)
load('EricModelDataSep14','GRpars','Dusp1FitParameters')
ModelGR.parameters(:,2) = num2cell([Dusp1FitParameters;GRpars']);

%% Associate Data with Different Instances of Model (1,10,100nm Dex)
GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
    '10','10',102,'GR Fit (10nM Dex)';...
    '100','100',103,'GR Fit (100nM Dex)'};

ModelGRparameterMap = cell(1,size(GRfitCases,1));
ModelGRfit = cell(1,size(GRfitCases,1));
for i=1:size(GRfitCases,1)
    ModelGRfit{i} = ModelGR.loadData('DUSP1_GR_dataframes/GR_ICC_3hr_Dex_total.csv',...
        {'nucGR','GRNorm'},...
        {'conc',GRfitCases{i,2}}); 
    ModelGRfit{i}.inputExpressions = {'IDex',[GRfitCases{i,2},'*exp(-gDex*t)*(t>0)']};
    ModelGRparameterMap(i) = {(1:7)};
end

%%  Combine all three GR models and fit with single parameter set.
fitOptions = optimset('Display','iter','MaxIter',1000);
combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap);
combinedGRModel = combinedGRModel.initializeStateSpaces;
GRpars = combinedGRModel.maximizeLikelihood(...
    GRpars, fitOptions);

%%  Make Plots of GR Fit Results
fignums = [111,121,GRfitCases{1,3},131;112,122,GRfitCases{2,3},132;113,123,GRfitCases{3,3},133];
combinedGRModel = combinedGRModel.updateModels(GRpars,true,fignums);
for i=1:size(GRfitCases,1)
    figure(GRfitCases{i,3}); 
    set(gca,'ylim',[0,18])
    title(GRfitCases{i,4})
    ylabel('Nuclear GR')
    xlabel('Time (min)')
    
    %  Update parameters in original models.
    ModelGRfit{i} = combinedGRModel.SSITModels{i};
end

%%  Extend Model to Include DUSP1 Activation, Production, and Degradation
ModelGRDusp = ModelGRfit{1};
ModelGRDusp.species = {'offGene';'onGene';'cytGR';'nucGR';'rna'};
ModelGRDusp.initialCondition = [2;0;24;1;5];
ModelGRDusp.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
    '(kcn0+kcn1*IDex)*cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR';...
    'kr*onGene';'gr*rna'};
ModelGRDusp.stoichiometry = [-1,1,0,0,0,0,0,0,0;...
                         1,-1,0,0,0,0,0,0,0;...
                         0,0,-1,1,1,-1,0,0,0;...
                         0,0,1,-1,0,0,-1,0,0;...
                         0,0,0,0,0,0,0,1,-1];
ModelGRDusp.useHybrid = true;
ModelGRDusp.hybridOptions.upstreamODEs = {'cytGR','nucGR'};
ModelGRDusp.solutionScheme = 'FSP';
ModelGRDusp.fspOptions.bounds = [0;0;0;2;2;400];
ModelGRDusp.fittingOptions.modelVarsToFit = 1:4;

%% Load and Associate with DUSP1 smFISH Data (100nM Dex Only)
Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
for i = 1:size(Dusp1FitCases,1)
    ModelDusp1Fit{i} = ModelGRDusp.loadData('DUSP1_GR_dataframes/DUSP1_3hr_Dex_100nM_total.csv',...
        {'rna','RNA_nuc'}); 
    ModelDusp1Fit{i}.inputExpressions = {'IDex',[Dusp1FitCases{i,2},'*exp(-gDex*t)*(t>0)']};    
    ModelDusp1parameterMap{i} = (1:4);
end
DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp.fittingOptions.modelVarsToFit,2}];

%%  Combine DUSP1 models and fit with single parameter set.
fitOptions = optimset('Display','iter','MaxIter',100);
fitOptions.suppressFSPExpansion = true; 
combinedDusp1Model = SSITMultiModel(ModelDusp1Fit,ModelDusp1parameterMap);
combinedDusp1Model = combinedDusp1Model.initializeStateSpaces;
DUSP1pars = combinedDusp1Model.maximizeLikelihood(...
    DUSP1pars, fitOptions);
ModelGRDusp.parameters(1:4,2) = num2cell(DUSP1pars);

%%  Make Plots of DUSP1 Fit Results
fignums = [211,221,201,231];
combinedDusp1Model = combinedDusp1Model.updateModels(DUSP1pars,true,fignums);
for i=1:size(Dusp1FitCases,1)
    figure(Dusp1FitCases{i,3}); 
    set(gca,'ylim',[0,150])
    title(Dusp1FitCases{i,4})
    ylabel('Nuclear DUSP1 mRNA')
    xlabel('Time (min)')
    
    %  Update parameters in original models.
    ModelDusp1Fit{i} = combinedDusp1Model.SSITModels{i};
end

%%  Predict DUSP1 Distributions under other Dex concentrations.
PredictionCases = {'10','10',301,'DUSP1 Prediction (10nM Dex)';...
                    '1','1',302,'DUSP1 Prediction (1nM Dex)';...
                    '0p3','0.3',303,'DUSP1 Prediction (0.3nM Dex)'};

fignums = [311,321,301,331;...
    312,322,302,332;...
    313,323,303,333];
ModelPred = cell(size(Dusp1FitCases,1),1);
for i=1:size(PredictionCases,1)
    ModelPred{i} = ModelGRDusp.loadData('DUSP1_GR_dataframes/DUSP1_3hr_Dex_TimeConcSweep_total.csv',...
        {'rna','RNA_nuc'},...
        {'conc',PredictionCases{i,2}}); 
    ModelPred{i}.inputExpressions = {'IDex',[PredictionCases{i,2},'*exp(-gDex*t)*(t>0)']};
    ModelPred{i}.makeFitPlot([],5,fignums(i,:))
    
    figure(fignums(i,3)); 
    set(gca,'ylim',[0,120])
    title(PredictionCases{i,4})
    ylabel('DUSP1 mRNA')
    xlabel('Time (min)')
end

%%  Predict DUSP1 Distributions under Dex Titration
DexConc = 10.^[-3,-2,-1,0,1,3,4];
DecConcStr = {'0.001','0.01','0.1','1','10','1000','10000'};

ModelPredDexTtr = cell(size(DexConc,1),1);
ModelPredDexTtrSoln = cell(size(DexConc,1),1);
for i=1:length(DexConc)
    ModelPredDexTtr{i} = ModelGRDusp.loadData('DUSP1_GR_dataframes/DUSP1_75min_ConcSweep_total.csv',...
        {'rna','RNA_nuc'},...
        {'conc',DecConcStr{i}}); 
    ModelPredDexTtr{i}.inputExpressions = {'IDex',[DecConcStr{i},'*exp(-gDex*t)*(t>0)']};
    ModelPredDexTtrSoln{i} = ModelPredDexTtr{i}.solve;

    DataHist = double(ModelPredDexTtr{i}.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor);
    PModel = double(ModelPredDexTtrSoln{i}.fsp{2}.p.sumOver([1,2]).data);

    MeanData(i) = DataHist*[0:length(DataHist)-1]'/sum(DataHist);
    MeanModel(i) = PModel'*[0:length(PModel)-1]'/sum(PModel);
    
    Mean2Data = DataHist*([0:length(DataHist)-1]').^2/sum(DataHist);
    Mean2Model = PModel'*([0:length(PModel)-1]').^2/sum(PModel);

    SigDat(i) = sqrt(Mean2Data-MeanData(i)^2);
    SigMod(i) = sqrt(Mean2Model-MeanModel(i)^2);
    
end
figure(401); clf; hold on
x = [DexConc,DexConc(end:-1:1)];
y = [MeanModel-SigMod,MeanModel(end:-1:1)+SigMod(end:-1:1)];
fill(x,y,[0.9,1,0.9])
plot(DexConc,MeanModel,'k','LineWidth',3); hold on;

errorbar(DexConc,MeanData,SigDat,'LineWidth',3,'LineStyle','none'); hold on;
plot(DexConc,MeanData,'s','MarkerSize',18,'MarkerFaceColor','b');

set(gca,'xscale','log','FontSize',16,'xlim',[9e-4,1.1e4])
title('Prediction of DUSP1 Expression')
xlabel('Dex Concentration at 75 min')
ylabel('DUSP1 Expression (nM)')

%%  Predict DUSP1 Distributions After Tryptolide
fignums = [511,521,501,531;...
    512,522,502,532;...
    513,523,503,533];

PredictionCases = {'20',502,'DUSP1 Prediction (t_{TPL} = 20 min)';...
                    '75',503,'DUSP1 Prediction (t_{TPL} = 75 min)';...
                    '180',504,'DUSP1 Prediction (t_{TPL} = 180 min)'};

ModelPredDexTpl = cell(length(tTPL),1);
ModelPredDexTplSoln = cell(length(tTPL),1);
for i=1:size(PredictionCases,1)
    ModelPredDexTpl{i} = ModelGRDusp.loadData('DUSP1_GR_dataframes/DUSP1_100nM_Dex_5uM_TPL_R1.csv',...
        {'rna','RNA_nuc'},...
        {'t_TPL',PredictionCases{i,1}}); 

    ModelPredDexTpl{i}.tSpan = sort(unique([ModelPredDexTpl{i}.tSpan,linspace(0,250,30)]));

    ModelPredDexTpl{i}.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
        '(kcn0+kcn1*IDex)*cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR';...
        'kr*onGene*Itrpt';'gr*rna'};

    ModelPredDexTpl{i}.inputExpressions = {'IDex','100*exp(-gDex*t)*(t>0)';...
        'Itrpt',['t<=',PredictionCases{i,1}]};
    
    ModelPredDexTpl{i}.makeFitPlot([],5,fignums(i,:))
    
    figure(fignums(i,3)); 
    set(gca,'ylim',[0,120])
    title(PredictionCases{i,3})
    ylabel('DUSP1 mRNA')
    xlabel('Time (min)')
end
%%




%%
% 
% if ~exist('ModelGRDuspPred10','var')
%     ModelGRDuspPred10 = ModelGRDusp;
%     ModelGRDuspPred10 = ModelGRDuspPred10.loadData('DUSP1_GR_dataframes/DUSP1_3hr_Dex_TimeConcSweep_total.csv',...
%         {'rna','RNA_nuc'},...
%         {'conc','10'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;
%     ModelGRDuspPred10.inputExpressions = {'IDex','10*exp(-gDex*t)*(t>0)'};
% end
% ModelGRDuspPred10.makeFitPlot([],5,[1,2,10,4])
% figure(10); set(gca,'ylim',[0,120])
% title('Dusp1 Prediction (10nM Dex)')
% ylabel('Dusp1 mRNA')
% xlabel('Time (min)')
% 
% if ~exist('ModelGRDuspPred1','var')
%     ModelGRDuspPred1 = ModelGRDusp;
%     ModelGRDuspPred1 = ModelGRDuspPred1.loadData('DUSP1_GR_dataframes/DUSP1_3hr_Dex_TimeConcSweep_total.csv',...
%         {'rna','RNA_nuc'},...
%         {'conc','1'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;
%     ModelGRDuspPred1.inputExpressions = {'IDex','1*exp(-gDex*t)*(t>0)'};
% end
% ModelGRDuspPred1.makeFitPlot([],5,[1,2,11,4])
% f = figure(11); set(gca,'ylim',[0,120])
% title('Dusp1 Prediction (1nM Dex)')
% ylabel('Dusp1 mRNA')
% xlabel('Time (min)')
% 
% if ~exist('ModelGRDuspPred0p3','var')
% ModelGRDuspPred0p3 = ModelGRDusp;
% ModelGRDuspPred0p3 = ModelGRDuspPred0p3.loadData('DUSP1_GR_dataframes/DUSP1_3hr_Dex_TimeConcSweep_total.csv',...
%     {'rna','RNA_nuc'},...
%     {'conc','0.3'}); % This would load the data assign x1 and x2 and condition on Rep_num = 1;
% ModelGRDuspPred0p3.inputExpressions = {'IDex','0.3*exp(-gDex*t)*(t>0)'};
% ModelGRDuspPred0p3.makeFitPlot([],5,[1,2,12,4])
% figure(12); set(gca,'ylim',[0,120])
% title('Dusp1 Prediction (0.3nM Dex)')
% ylabel('Dusp1 mRNA')
% xlabel('Time (min)')
% 

%% Create and estimate parameters for combined models
%%      Example 0 -- single model.
% This is a simple example, where we only fit one model to a single data
% set.
singleModel = SSITMultiModel({Model},{[1:7]});
singleModel = singleModel.initializeStateSpaces;
allParsSingle = ([Model.parameters{:,2}]);
allParsSingle = singleModel.maximizeLikelihood(...
    allParsSingle, fitOptions, fitAlgorithm);
singleModel = singleModel.updateModels(allParsSingle);

% Copy the parameters back into Model1 and Model2 so we can reuse them later.
Model.parameters = singleModel.SSITModels{1}.parameters;
ModelGR.parameters = singleModel.SSITModels{1}.parameters;

%%      Example 1 -- Adding new Model+Data to an existing multimodel
% This is how one adds a second model/data combination.  In this case the
% parameters of the new model are completely independent of the parameter
% set for the first model.
combinedModel = singleModel.addModel({ModelGR},{[8:14]});
combinedModel = combinedModel.initializeStateSpaces;
allParsCombined = ([Model.parameters{:,2},[ModelGR.parameters{:,2}]]);
allParsCombined = combinedModel.maximizeLikelihood(...
    allParsCombined, fitOptions, fitAlgorithm);
combinedModel = combinedModel.updateModels(allParsCombined);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting a single model independently, then it is more efficient to fit
% it separately.

%%      Example 2 -- completely independent parameters.
% Here is how we can create the combined model in one shot.
combinedModelIndependent = SSITMultiModel({Model,ModelGR},{[1:7],[8:14]});
combinedModelIndependent = combinedModelIndependent.initializeStateSpaces;
allParsIndepdendent = ([Model.parameters{:,2},[ModelGR.parameters{:,2}]]);
allParsIndepdendent = combinedModelIndependent.maximizeLikelihood(...
    allParsIndepdendent, fitOptions, fitAlgorithm);
combinedModelIndependent = combinedModelIndependent.updateModels(allParsIndepdendent);

%%      Example 3 -- completely dependent parameters.
% Here is an example of how a single set of parameters can be used for both
% models and data sets. In the following we make a joint model where both
% Model1 and Model2 use the parameters [1:7].
combinedModelDependent = SSITMultiModel({Model,ModelGR},{[1:7],[1:7]});
combinedModelDependent = combinedModelDependent.initializeStateSpaces;
allParsDependent = ([Model.parameters{:,2}]);
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);
combinedModelDependent = combinedModelDependent.updateModels(allParsDependent);

% Note: This example is shown for illustration purposes only.  Usually, if 
% one is fitting two replicas of the exact same experiment, then it is
% more efficient to combine the data from both replicas and fit them at the
% same time, e.g. to combined all replicas into one set, simply loas the
% data as follows:
% Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});

%%      Example 4 -- mixed parameters.
% Sometimes it is desirable to only let some parameters change from
% condition to condition.  In this example both Model1 and Model2 use the
% same parameters [1-4], but parameters [5:7] are only for Model1 and
% [8:10] are only for Model2.
combinedModelMixed= SSITMultiModel({Model,ModelGR},{[1:7],[1:4,8:10]});
combinedModelMixed = combinedModelMixed.initializeStateSpaces;
allParsMixed = ([Model.parameters{:,2},ModelGR.parameters{5:7,2}]);
allParsMixed = combinedModelMixed.maximizeLikelihood(...
    allParsMixed, fitOptions, fitAlgorithm);

%%      Example 5 -- constrained parameters.
% it is also often helpful to place constraints on parameters since it can
% be expectd that cartain parameters should not change that much from one
% experiment to another.
constraint = @(x)sum((x(5:7)-x(8:10)).^2);
combinedModelConstrained = SSITMultiModel({Model,ModelGR},{[1:7],[1:4,8:10]},constraint);
combinedModelConstrained = combinedModelConstrained.initializeStateSpaces;
allParsConstrained = ([Model.parameters{:,2},ModelGR.parameters{5:7,2}]);
allParsConstrained = combinedModelConstrained.maximizeLikelihood(...
    allParsConstrained, fitOptions, fitAlgorithm);