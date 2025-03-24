%% Additional Codes for FIM and PDO analyses (ALEX Paper)
%%  STEP PDO -- Explore Optimal Designs Using FIM Add FIM to MH UQ Plots 
%%    STEP PDO.A. -- FIM Analyses
%%      STEP PDO.A.1. -- Plot UQ from FIM compared to MH
close all
clear
addpath(genpath('../../src'));
addpath(genpath('tmpPropensityFunctions'));

loadPrevious = true;
%savedWorkspace = 'workspaceJuly24';
savedWorkspace = 'workspaceOct22_2024';
%savedWorkspace = 'workspaceDec9_2025';

if loadPrevious
    load(savedWorkspace);
end

figNew = figure;

GRDusp1_log10PriorStd = log10PriorStd(1:13);

fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
    diag(GRDusp1_log10PriorStd.^2));
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

%%      STEP PDO.A.2. -- Find optimal experiment design (same number of cells)
nTotal = sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOpt = ModelGRDusp100nM.optimizeCellCounts(fimResults,nTotal,'tr[1:4]'); % min. inv determinant <x^{-1}> (all other parameters are known and fixed)
nCellsOptAvail = min(nCellsOpt,ModelGRDusp100nM.dataSet.nCells')
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));
fimOptAvail = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptAvail,diag(GRDusp1_log10PriorStd.^2));
figOpt = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimTotal],'log',[],figOpt);
figOptAvail = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOptAvail,fimTotal],'log',[],figOptAvail);
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


%%    STEP PDO.B. -- PDO Calculations
%%      STEP PDO.B.1. -- Calibrate PDO from Multi-Modal Experimental Data
% Calibration the PDO from empirical data. Here, the number of spots has
% been measured using different assays in data columns 'nTotal' for the
% 'true' data set and in the columns 'nSpots0' for a different label or
% 'intens1' for the integrated intensity.  We calibrate two different PDOs
% for this case. In both cases, we assume an 'AffinePoiss' PDO where the
% obervation probability is a Poisson distribution where the mean value is
% affine linearly related to the true value: P(y|x) = Poiss(a0 + a1*x);
ModelPDOSpots = ModelGRDusp100nM.calibratePDO('../../ExampleData/pdoCalibrationData.csv',...
    {'rna'},{'nTotal'},{'nSpots0'},'AffinePoiss',true);

%% STEP PDO.B.1.a -- Calibrate PDO from single (and 5-9) z-stacks
ModelPDO_single_Z = ModelGRDusp100nM;
ModelPDO_5thru9_Zs = ModelGRDusp100nM;
ModelPDO_single_Z = ModelPDO_single_Z.calibratePDO('../AlexP/SpotCountingByCellAndByZ/merged_dataframe_A549_DUSP1_100nM_10min_062723.csv',...
    {'rna'},{'Spot_Count'},{'Z_7'},'AffinePoiss',true);
ModelPDO_5thru9_Zs = ModelPDO_5thru9_Zs.calibratePDO('../AlexP/SpotCountingByCellAndByZ/merged_dataframe_A549_DUSP1_100nM_10min_062723.csv',...
    {'rna'},{'Spot_Count'},{'Z_5_9'},'AffinePoiss',true);

%%      STEP PDO.B.2. -- Calibrate PDO from Eric's DUSP1 Intensity Data.
ModelPDOIntensEric = ModelGRDusp100nM;
ModelPDOIntensEric = ModelPDOIntensEric.calibratePDO('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
    {'rna'},{'RNA_DUSP1_nuc'},{'Nuc_DUSP1_avg_int_tot'},'AffinePoiss',true,[1,230,0.5]);

%%    STEP PDO.C. -- FIM + PDO Analyses
%%      STEP PDO.C.1. -- Analyze FIM with PDO for MCP/smFISH
 fimsPDOSpot = ModelPDOSpots.computeFIM([],'log');
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));

nCellsOptPDOspots = ModelPDOSpots.optimizeCellCounts(fimsPDOSpot,nTotal,'tr[1:4]');


figMCP = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimPDOSpots,fimTotal,fimOpt],'log',[],figMCP);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[0,0,0];   % MH - black
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];   % MLE - black
        CH(2).LineWidth = 3;
        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
        CH(3).LineWidth = 3;
        CH(4).Color=[0,1,1];   % fimTotal - blue 
        CH(4).LineWidth = 3;
        CH(5).Color=[1,0,1];   % fimOpt - magenta
        CH(5).LineWidth = 3;
    end
end

%%      STEP PDO.C.1.a. -- Analyze FIM with PDO for single Z and Zs 5-9
fimsPDO_single_Z = ModelPDO_single_Z.computeFIM([],'log');
fimPDO_single_Z = ModelPDO_single_Z.evaluateExperiment(fimsPDO_single_Z,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));

nCellsOptPDO_single_Zs = ModelPDO_single_Z.optimizeCellCounts(fimsPDO_single_Z,nTotal,'tr[1:4]');

fig_single_Z = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimPDO_single_Z,fimTotal,fimOpt],'log',[],fig_single_Z);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[0,0,0];   % MH - black
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];   % MLE - black
        CH(2).LineWidth = 3;
        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
       CH(3).LineWidth = 3;
        CH(4).Color=[0,1,1];   % fimTotal - blue 
        CH(4).LineWidth = 3;
        CH(5).Color=[1,0,1];   % fimOpt - magenta
        CH(5).LineWidth = 3;
    end
end

fimsPDO_5thru9_Zs = ModelPDO_5thru9_Zs.computeFIM([],'log');
fimPDO_5thru9_Zs = ModelPDO_5thru9_Zs.evaluateExperiment(fimsPDO_5thru9_Zs,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));

nCellsOptPDO_5thru9_Zs = ModelPDO_5thru9_Zs.optimizeCellCounts(fimsPDO_5thru9_Zs,nTotal,'tr[1:4]');

fig_Zs_5thru9 = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimPDO_5thru9_Zs,fimTotal,fimOpt],'log',[],fig_Zs_5thru9);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[0,0,0];   % MH - black
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];   % MLE - black
        CH(2).LineWidth = 3;
        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
        CH(3).LineWidth = 3;
        CH(4).Color=[0,1,1];   % fimTotal - blue 
        CH(4).LineWidth = 3;
        CH(5).Color=[1,0,1];   % fimOpt - magenta
        CH(5).LineWidth = 3;
    end
end

%%      STEP PDO.C.2. -- Analyze FIM with PDO for Intensity only
fimsPDOIntens = ModelPDOIntensEric.computeFIM([],'log');
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));

nCellsOptPDOintens = ModelPDOSpots.optimizeCellCounts(fimsPDOIntens,nTotal,'tr[1:4]');
 
fimPDOIntensAvail = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptAvail,diag(GRDusp1_log10PriorStd.^2));

figInt = figure; clf;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimPDOSpots,fimPDOIntens,fimPDOIntensAvail],'log',[],figInt);
% figNew = figure; clf;
% ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimTotal,fimPDOIntensAvail],'log',[],figNew);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[0,0,0];   % MH - black
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];   % MLE - black
        CH(2).LineWidth = 3;
        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
        CH(3).LineWidth = 3;
        CH(4).Color=[0,1,1];   % fimPDOIntens - blue
        CH(4).LineWidth = 3;
        CH(5).Color=[1,0,1];   % fimPDOIntensAvail - magenta
        CH(5).LineWidth = 3;
    end
end

% figNew = figure; clf;
% ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimOpt,fimTotal,fimPDOIntens],'log',[],figNew);
%for i = 1:3
%    for j = i:3
%        subplot(3,3,(i-1)*3+j)
%        CH = get(gca,'Children');
%        CH(1).Color=[0,0,0];   % MH - black
%        CH(1).LineWidth = 3;
%        CH(2).Color=[0,0,0];   % MLE - black
%        CH(2).LineWidth = 3;
%        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
%        CH(3).LineWidth = 3;
%        CH(4).Color=[0,1,1];   % fimTotal - blue
%        CH(4).LineWidth = 3;
%        CH(5).Color=[1,0,1];   % fimPDOIntens - magenta
%        CH(5).LineWidth = 3;
%    end
%end

%%      STEP PDO.C.3. -- Analyze FIM with PDO for Intensity only more cells
nTimes = 3.71;
fimPDOIntens2x = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt*nTimes,diag(GRDusp1_log10PriorStd.^2));
det(fimOpt{1}(1:4,1:4))/det(fimPDOIntens2x{1}(1:4,1:4))

figIntMoreCells = figure;
ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[fimPDOSpots,fimTotal,fimPDOIntens2x],'log',[],figIntMoreCells);

for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[0,0,0];   % MH - black
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];   % MLE - black
        CH(2).LineWidth = 3;
        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
        CH(3).LineWidth = 3;
        CH(4).Color=[0,1,1];   % fimTotal - blue
        CH(4).LineWidth = 3;
        CH(5).Color=[1,0,1];   % fimPDOIntens2x - magenta
        CH(5).LineWidth = 3;
    end
end


%%    STEP PDO.D. -- Validation of FIM Predictions
%%      STEP PDO.D.1. -- Fit to FIM selected time points
ModelGRDusp100nM_FIMDesign = ModelGRDusp100nM;
% Set the fitting routine to only consider the time points selected by the
% FIM analysis:
ModelGRDusp100nM_FIMDesign.fittingOptions.timesToFit = nCellsOpt>0;
% Refit the model, but now with only those time points.
DUSP1parsFIMDesign = DUSP1pars;
fitOptions = optimset('Display','iter','MaxIter',200);
for i = 1:1
    DUSP1parsFIMDesign = ModelGRDusp100nM_FIMDesign.maximizeLikelihood(...
        DUSP1parsFIMDesign, fitOptions);
    ModelGRDusp100nM_FIMDesign.parameters(1:4,2) = num2cell(DUSP1parsFIMDesign);
    save('EricModelDusp1_MMDex_PDO','GRpars','DUSP1pars','DUSP1parsFIMDesign') 
end

%%      STEP PDO.D.2. -- Plot fit and predicitons using FIM suggested conditions.
showCases = [1,1,1,1];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsFIMDesign,Dusp1FitCases,showCases)

%%      STEP PDO.D.3. -- Fit the DISTORTED intensity measurements at ALL times.
% Here, I will load the intensity data for the nuclear DUSP1.  then I will
% attempt to identify the model from just that data and at just the times
% selected by the FIM.
ModelPDOIntensEric.dataSet = [];
ModelPDOIntensEric = ModelPDOIntensEric.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','Nuc_DUSP1_avg_int_tot'},...
        {'Dex_Conc','100'}); 
load('EricModelDusp1_MMDex_PDO','DUSP1pars') 

ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1pars);

fitIters=3;
for i = 1:fitIters
    DUSP1parsIntensity = ModelPDOIntensEric.maximizeLikelihood(...
        DUSP1pars, fitOptions);
    ModelPDOIntensEric.parameters(1:4,2) = num2cell(DUSP1parsIntensity);
    save('EricModelDusp1_MMDex_PDO','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsIntensity') 
end
%%        STEP PDO.D.3.a. -- Plot predictions when fit to distorted data at ALL times.
showCases = [1,1,1,1];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsIntensity,Dusp1FitCases,showCases)

%%      STEP PDO.D.4. -- Fit the DISTORTED intensity measurements at FIM selected times.
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
load('EricModelDusp1_MMDex_PDO','DUSP1parsFIMDesign') 
ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesign);

for i = 1:5
    DUSP1parsFIMDesignIntensity = ModelPDOIntensEricFIM.maximizeLikelihood(...
        DUSP1parsFIMDesign, fitOptions);
    ModelPDOIntensEricFIM.parameters(1:4,2) = num2cell(DUSP1parsFIMDesignIntensity);
    save('EricModelDusp1_MMDex_PDO','GRpars','DUSP1pars','DUSP1parsFIMDesign','DUSP1parsFIMDesignIntensity') 
end
%%        STEP PDO.D.4.a. -- Plot predictions when fit to distorted data at FIM times.
showCases = [1,1,1,1];
makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1parsFIMDesignIntensity,Dusp1FitCases,showCases)

%%    STEP PDO.E. -- Plot Information vs. Expt Design, PDO, and Number of Cells
fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,diag(GRDusp1_log10PriorStd.^2));
fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));
fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOpt,diag(GRDusp1_log10PriorStd.^2));

barWithOriginalNumber = [det(fimTotal{1}(1:4,1:4)),det(fimOpt{1}(1:4,1:4)),det(fimPDOSpots{1}(1:4,1:4)),det(fimPDOIntens{1}(1:4,1:4))];

nCellVec = logspace(3,5,20);
nCellsOrigRat = ModelGRDusp100nM.dataSet.nCells/sum(ModelGRDusp100nM.dataSet.nCells);
nCellsOptRat = nCellsOpt/sum(nCellsOpt);

cols = {'k','c','b','r'};  % fimOrig (black), fimOpt (cyan), fimPDOSpots (blue), fimPDOIntens (red) 
fimDetVsNumberMAt=[];
for i = 1:length(nCellVec)
    fimOrig = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOrigRat*nCellVec(i),diag(GRDusp1_log10PriorStd.^2));
    fimOpt = ModelGRDusp100nM.evaluateExperiment(fimResults,nCellsOptRat*nCellVec(i),diag(GRDusp1_log10PriorStd.^2));
    fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOptRat*nCellVec(i),diag(GRDusp1_log10PriorStd.^2));
    fimPDOIntens = ModelPDOIntensEric.evaluateExperiment(fimsPDOIntens,nCellsOptRat*nCellVec(i),diag(GRDusp1_log10PriorStd.^2));
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
