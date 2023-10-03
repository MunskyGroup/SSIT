%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
clear all
close all
clc
addpath('../CommandLine');
addpath('../EricModel/DUSP1_GR_dataframes');
addpath(genpath('../src'));
addpath('../tensor_toolbox-v3.2.1')
%% Model Analysis Options
% Use previous saved results
useSavedModel = true;
useSavedModelMH = false;
useSavedModelSensitivty = true;
useSavedModelMLE = false; 

useManualModelCreation = true; % Create model from scratch
fitModel = true; % Fit model with fminsearch
runMetHast = true; % Calculation of Metropolis Hastings 

calcMLE = false; % MLE calculation for FIM verification
runOptTriptolideExp = false; % optimal triptolide application calc

showFigures = true;

% Save Model, MHresults, sens, and MLE results
% Previous results will not be overwritten
saveModelOutput = true; 
saveMHOutput = true; 
saveSensOutput = true; 
saveMLEOutput = true; 

%% EGRNT Dusp1 model
if useManualModelCreation
    Model = SSIT;
    Model.species = {'x1';'x2';'x3'};  % x1: nuclear GR x2: number of on alleles,  x3: mRNA 
    Model.initialCondition = [0;0;0];
    Model.propensityFunctions = {'(kcn0+kcn1*IDex)';'knc*x1';...
                                 'kon*x1*(2-x2)';'koff*x2';'kr*x2';'gr*x3'};
    
    stoich = [1,0,0;   % GR enters nucleus
             -1,0,0;   % GR exist nucleus / degrades
              0,1,0;   % Gene allele switching to on state
              0,-1,0;  % Gene allele switching to off state
              0,0,1;   % mRNA production
              0,0,-1]; % mRNA degredation
    
    Model.stoichiometry = transpose(stoich);
    
    % Input expression for time varying signals
    Model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'}; % Expression for Dex stimulation
    
    % Defining the values of each parameter used
    Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
                         'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
    
    Model.fspOptions.initApproxSS = true;
end

%% Load saved mHast and Sensitivity
% Only run if the result mat files are already created
if useSavedModel
%     Model = load('EGRNT_model.mat').Model;
    Model = load('complex_dusp1_model.mat').Model;

end

if useSavedModelMH
    mhResults = load('EGRNT_mhast.mat').mhResults;
%     mhResults = load('complex_dusp1_mhast.mat').mhResults;

end

if useSavedModelSensitivty
%     sensSoln = load('EGRNT_sens.mat').sensSoln;
    sensSoln = load('complex_dusp1_sens.mat').sensSoln;

end

%% Solve the model using the FSP
if fitModel
    Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'}); % Load experimental data set
    Model.fittingOptions.modelVarsToFit = 1:8; % Picking parameters 1-8 for fitting
    
    for i=1:5
      Model.solutionScheme = 'FSP';
      Model.fspOptions.fspTol = 1e-4;
      Model.fspOptions.verbose = 0;
      Model.fspOptions.bounds=[];
      [fspSoln,Model.fspOptions.bounds] = Model.solve;
      Model.fspOptions.bounds
    
      % Load and Fit smFISH Data
      Model.fspOptions.fspTol = inf;
      Model.fittingOptions.modelVarsToFit = 1:8;
      fitOptions = optimset('Display','iter','MaxIter',500);
      Model.parameters(1:8,2) = num2cell(Model.maximizeLikelihood([],fitOptions));
    end
end
% figure slide 5 and 6
if showFigures
    Model.makeFitPlot;
end
%% Metropolis Hastings to Quantify Parameter Uncertainty
if runMetHast
        
    Model.fittingOptions.modelVarsToFit = 1:4; % Picking parameters 1-4 for fitting
    
    % MH options set to save 15000 samples, Fisrt 100 samples discarded, and
    % every other sample taken being saved
    MHOptions = struct('numberOfSamples',15000,'burnin',100,'thin',1,...
    'useFIMforMetHast',true,'suppressFSPExpansion',true);
    
    % The parameter set that maximizes the likelihood function is saved with
    % the MH samples
    [bestParsFound,~,mhResults] = Model.maximizeLikelihood([Model.parameters{Model.fittingOptions.modelVarsToFit,2}]',...
    MHOptions,'MetropolisHastings');
    Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
end

Model.fittingOptions.modelVarsToFit = 1:8;

% figure slide 9
if showFigures
    Model.plotMHResults(mhResults);
end

%% Calculate CME Sensitivity to Parameter Variations
Model.sensOptions.solutionMethod = 'finiteDifference';
Model.solutionScheme = 'fspSens';
Model.fspOptions.fspTol = 1e-6;
[sensSoln,Model.fspOptions.bounds] = Model.solve;
%% Maximum likelihood Estimator Verification
% Simmulate data with known parameters and re-fit simmulated data to get MLE 
starttime=-1500;
EGRNTssa=Model;
EGRNTssa.initialTime = starttime;
EGRNTssa.tSpan = [starttime,EGRNT_Model_check.tSpan];
EGRNTssa.ssaOptions.useTimeVar=true;
EGRNTssa.ssaOptions.signalUpdateRate=1;
EGRNTssa.solutionScheme = 'SSA';  % Set solution scheme to SSA.
EGRNTssa.ssaOptions.Nexp = 200; 
EGRNTssa.ssaOptions.nSimsPerExpt = 200;
EGRNTssa.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
EGRNTssa.solve([],'EGRNT_SSA.csv'); 

%% Find MLE for simulated data set.

EGRNT_Model_check.fittingOptions.modelVarsToFit = [2,3];
%nFrePars = length(comp_Model.fittingOptions.modelVarsToFit)
nFrePars=2;
EGRNT_MLE = zeros(1,nFrePars,EGRNTssa.ssaOptions.Nexp);
fMLE = inf(1,nFrePars,EGRNTssa.ssaOptions.Nexp);
B = EGRNT_Model_check.loadData('EGRNT_SSA.csv',{'x1','exp1_s1';'x2','exp1_s2';'x3','exp1_s3'});

for iExp = 1:EGRNTssa.ssaOptions.Nexp
    
    % GR & Dusp1  observed
    B = EGRNT_Model_check.loadData('EGRNT_SSA.csv',{'x2',['exp',num2str(iExp),'_s2'];'x3',['exp',num2str(iExp),'_s3']}); % Link non-distorted data.

    fitOptions = optimset;
    fitOptions.MaxIter = 500;
    fitOptions.Display = 'iter';
    
    B.fittingOptions.timesToFit = [false,ones(1,length(EGRNT_Model_check.tSpan),'logical')];
    B.tSpan = B.tSpan(2:end);
%     B.fittingOptions.modelVarsToFit = [2,3];
     if iExp==1
         x0 = [B.parameters{B.fittingOptions.modelVarsToFit,2}]';
     else
         x0 = squeeze(EGRNT_MLE(1,:,iExp-1)); 
     end
    [EGRNT_MLE(1,:,iExp),fMLE(1,:,iExp)] = B.maximizeLikelihood(x0,fitOptions);
end

% Solve for Fisher Information Matrix at all Time Points
B.fittingOptions.modelVarsToFit = [1:4];

B.solutionScheme = 'FSP';
B.fspOptions.fspTol = 1e-6;
B.fspOptions.bounds=[];
[fspSoln,B.fspOptions.bounds] = B.solve;

B.fspOptions.fspTol = inf;
B.solutionScheme = 'fspSens';
sensSoln = B.solve(fspSoln.stateSpace);

B.pdoOptions.unobservedSpecies = {'x1','x2'};
fims_B = B.computeFIM(sensSoln.sens);
FIM_B = B.evaluateExperiment(fims_B,B.dataSet.nCells);


% Make Plots (figure slides 14 and 15)
fimFreePars = FIM_B(B.fittingOptions.modelVarsToFit,B.fittingOptions.modelVarsToFit);
if showFigures
    B.makeMleFimPlot(squeeze(EGRNT_MLE(1,:,:)),fimFreePars,[4,3],0.95)
end

%% Experiment Design Comparison
% Computation of each the FIM of each experiment design

%% Using all cells to measure DUSP1 expression (smFISH)
Model.pdoOptions.unobservedSpecies = {'x1','x2'}; % Gene state and nuclear GR are unobserved
fims_dusp = Model.computeFIM(sensSoln.sens); 
FIM_dusp = Model.evaluateExperiment(fims_dusp,Model.dataSet.nCells);

nTotal_dusp = sum(Model.dataSet.nCells);
nCellsOpt_dusp = Model.optimizeCellCounts(fims_dusp,nTotal_dusp,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_dusp_case2 = Model.optimizeCellCounts(fims_dusp,nTotal_dusp,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt_dusp = Model.evaluateExperiment(fims_dusp,nCellsOpt_dusp);
fimOpt_dusp_case2 = Model.evaluateExperiment(fims_dusp,nCellsOpt_dusp_case2);


%% Measuring Nuclear Glucocotoroid (ICC)
Model.pdoOptions.unobservedSpecies = {'x3','x2'}; % Gene state and mRNA are unobserved
fims_GR = Model.computeFIM(sensSoln.sens);
FIM_GR = Model.evaluateExperiment(fims_GR,Model.dataSet.nCells);

nTotal_GR = sum(Model.dataSet.nCells);
nCellsOpt_GR = Model.optimizeCellCounts(fims_GR,nTotal_GR,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_GR_case2 = Model.optimizeCellCounts(fims_GR,nTotal_GR,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt_GR = Model.evaluateExperiment(fims_GR,nCellsOpt_GR);
fimOpt_GR_case2 = Model.evaluateExperiment(fims_GR,nCellsOpt_GR_case2);

%% Seperate experiments to measure DUSP1 and GR (separate smFISH and ICC experiments)
% The same population of cells as measuring just mRNA is used. Half of the
% cells are measuring nuclear GR while the other half measures mRNA
fims_dusp_and_gr_separate = [fims_dusp;fims_GR];

% This will fix the number of cells at half the original number for each
% experiment.
FIM = (FIM_GR+FIM_dusp)/2;

% Here we will also allow twice as many cells to be measured with this
% experiment in order to match the case above.
nTotal = sum(Model.dataSet.nCells);
nCellsOpt_sep = Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_sep_case2 = Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt = Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep);
fimOpt_case2 = Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep_case2);

%% Simultaneous measurement of DUSP1 and GR (smFISH and ICC experiments in the same population of cells)
Model.pdoOptions.unobservedSpecies = {'x2'}; % Gene state is unobserved
fims_both = Model.computeFIM(sensSoln.sens);
FIM_both = Model.evaluateExperiment(fims_both,Model.dataSet.nCells);

nTotal_both = sum(Model.dataSet.nCells);
nCellsOpt_both = Model.optimizeCellCounts(fims_both,nTotal_both,'TR[1:4]');% In Case 1, use "TR[1:4]" command
nCellsOpt_both_case2 = Model.optimizeCellCounts(fims_both,nTotal_both,'[1:4]');% In Case 2, use "[1:4]" command

fimOpt_both = Model.evaluateExperiment(fims_both,nCellsOpt_both);
fimOpt_both_case2 = Model.evaluateExperiment(fims_both,nCellsOpt_both_case2);

%% Plot for Optimal Cell allocation measurement for each experiment design
timepoints = {'0min','10min','20min','30min','40min','50min','60min','75min','90min','120min','150min','180min'};

% case1 figure slide 16
figure()
designs = [Model.dataSet.nCells;nCellsOpt_dusp;nCellsOpt_sep(1:12);nCellsOpt_both];
designs_b = [Model.dataSet.nCells;nCellsOpt_dusp;nCellsOpt_sep(1:12)+nCellsOpt_sep(13:24);nCellsOpt_both];


g = bar(designs_b','barWidth',1);
hold on
h = bar(designs','barWidth',1);
hold off

set(gca,'XTickLabel',timepoints,'FontSize',20)
xlabel('Measurement timepoints')
ylabel('Number of Cell')
title('Case 1 (GR parameters known)')
legend({'','',append('Seperate GR Experiments Optimized (',int2str(sum(nCellsOpt_sep(13:24))/sum(nCellsOpt_sep)*100),'% of cells)')...
  ,'','Equal cells at each timepoint','DUSP1 Optimized',append('Seperate DUSP1 Experiments Optimized (',int2str(sum(nCellsOpt_sep(1:12))/sum(nCellsOpt_sep)*100),'% of cells)'),'Simutaneous DUSP1/GR Experiments Optimized',''},'FontSize',20)

% case2 figure slide 17
figure()
designs = [Model.dataSet.nCells;nCellsOpt_dusp_case2;nCellsOpt_sep_case2(1:12);nCellsOpt_both_case2];
designs_b = [Model.dataSet.nCells;nCellsOpt_dusp_case2;nCellsOpt_sep_case2(1:12)+nCellsOpt_sep_case2(13:24);nCellsOpt_both_case2];

g = bar(designs_b','barWidth',1);
hold on
h = bar(designs','barWidth',1);
hold off

set(gca,'XTickLabel',timepoints,'FontSize',20)
xlabel('Measurement timepoints')
ylabel('Number of Cell')
title('Case 2 (GR parameters not known)')
legend({'','',append('Seperate GR Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(13:24))/sum(nCellsOpt_sep_case2)*100),'% of cells)')...
  ,'','Equal cells at each timepoint','DUSP1 Optimized',append('Seperate DUSP1 Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(1:12))/sum(nCellsOpt_sep_case2)*100),'% of cells)'),'Simutaneous DUSP1/GR Experiments Optimized',''},'FontSize',20)

%% Comparing FIM of each experiment design (Case 1)
% figure slide 18
y = [det(inv(FIM_dusp(1:4,1:4))),det(inv(fimOpt_dusp(1:4,1:4)));...
  det(inv(FIM(1:4,1:4))),det(inv(fimOpt(1:4,1:4))); det(inv(FIM_both(1:4,1:4))), det(inv(fimOpt_both(1:4,1:4)))];
x = categorical({'Dusp1', 'Seperate GR and Dusp1 experiment', 'Simultaneous GR and Dusp1 experiment'});
bar(x,y,'barWidth',1)
set(gca, 'YScale', 'log','ylim',[1e-36,1e-15])
ax=gca;
ax.FontSize = 20;
legend({'FIM','FIM optimized'},'FontSize', 24)
ylabel('|COV|')
title('Case 1 (GR parameters known)')

%% Comparing FIM of each experiment design (Case 2)
% figure slide 19
ev = [eye(4),zeros(4)];
y = [det(ev*inv(FIM_dusp(1:8,1:8))*ev'),det(ev*inv(fimOpt_dusp_case2(1:8,1:8))*ev');...
  det(ev*inv(FIM(1:8,1:8))*ev'),det(ev*inv(fimOpt_case2(1:8,1:8))*ev'); det(ev*inv(FIM_both(1:8,1:8))*ev'), det(ev*inv(fimOpt_both_case2(1:8,1:8))*ev')];
x = categorical({'Dusp1', 'Seperate GR and Dusp1 experiment', 'Simultaneous GR and Dusp1 experiment'});
bar(x,y,'barWidth',1)
set(gca, 'YScale', 'log','ylim',[1e-25,1e-8])
ax2=gca;
ax2.FontSize = 25;
legend({'FIM','FIM optimized'},'FontSize', 24)
ylabel('|COV|')
title('Case 2 (GR parameters not known')
nTotal_all = sum(Model.dataSet.nCells);

%% Effect of additional cell measurements on information gained
  n = logspace(0,8,1000);
  for i=1:length(n)
     inv_FIM_case1_dusp(i) = det(inv((n(i)/nTotal_all)*FIM_dusp(1:4,1:4)));
     inv_FIM_case1_dusp_opt(i) = det(inv((n(i)/nTotal_all)*fimOpt_dusp(1:4,1:4)));

     inv_FIM_case1_sep(i) = det(inv((n(i)/nTotal_all)*FIM(1:4,1:4)));
     inv_FIM_case1_sep_opt(i) = det(inv((n(i)/nTotal_all)*fimOpt(1:4,1:4)));

     inv_FIM_case1_both(i) = det(inv((n(i)/nTotal_all)*FIM_both(1:4,1:4)));
     inv_FIM_case1_both_opt(i) = det(inv((n(i)/nTotal_all)*fimOpt_both(1:4,1:4)));

     inv_FIM_case2_dusp(i) = det(ev*inv((n(i)/nTotal_all)*FIM_dusp(1:8,1:8))*ev');
     inv_FIM_case2_dusp_opt(i) = det(ev*inv((n(i)/nTotal_all)*fimOpt_dusp_case2(1:8,1:8))*ev');

     inv_FIM_case2_sep(i) = det(ev*inv((n(i)/nTotal_all)*FIM(1:8,1:8))*ev');
     inv_FIM_case2_sep_opt(i) = det(ev*inv((n(i)/nTotal_all)*fimOpt_case2(1:8,1:8))*ev');

     inv_FIM_case2_both(i) = det(ev*inv((n(i)/nTotal_all)*FIM_both(1:8,1:8))*ev');
     inv_FIM_case2_both_opt(i) = det(ev*inv((n(i)/nTotal_all)*fimOpt_both_case2(1:8,1:8))*ev');

  end
  lnwdt=2;
    
  % figure slide 22
  if showFigures
      figure()% Case 1 with optimization
      plot(n,inv_FIM_case1_dusp,'b',n,inv_FIM_case1_sep,'g',n,inv_FIM_case1_both,'m',...
          n,inv_FIM_case1_dusp_opt,'--b',n,inv_FIM_case1_sep_opt,'--g',n,inv_FIM_case1_both_opt,'--m','LineWidth',lnwdt)
      yline(9.97e-25,'k','LineWidth',lnwdt);
      set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-35,1e-10],'xlim',[0,1e8],'FontSize', 24)
      xlabel('Number of Cells Measured')
      ylabel('|COV|')
      legend({'DUSP1 experiment','Seperate GR and DUSP1 experiment','Simultaneous GR and DUSP1 experiment','','','','Desired Information for Experiments'},'FontSize',20);
      title('Case 1')
  end
  
  if showFigures
      % figure slide 23
      figure()%Case 2 with optimization
      plot(n,inv_FIM_case2_dusp,'b',n,inv_FIM_case2_sep,'g',n,inv_FIM_case2_both,'m',...
          n,inv_FIM_case2_dusp_opt,'--b',n,inv_FIM_case2_sep_opt,'--g',n,inv_FIM_case2_both_opt,'--m','LineWidth',lnwdt)
      yline(1e-14,'k','LineWidth',lnwdt);
      set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-20,1e-10],'xlim',[0,1e8],'FontSize', 24)
      xlabel('Number of Cells Measured')
      ylabel('|COV|')
      legend({'DUSP1 experiment','Seperate GR and DUSP1 experiment','Simultaneous GR and DUSP1 experiment','','','','Desired Information for Experiments'},'FontSize',20);
      title('Case 2')
  end


%% Triptolide Experiment
if runOptTriptolideExp
    ModelTript = EGRNT_Model_check;
    ModelTript.propensityFunctions(5) = {'kr*x2*Itrypt'};
    ModelTript.inputExpressions(2,:) = {'Itrypt','(t<tpt)'};
    
    
    ModelTript.sensOptions.solutionMethod = 'finiteDifference';
    
    tpt_array = [0:5:40,42:2:180]; % Range of times to calculate FIM
    sensSoln = cell(1,length(tpt_array));
    for iTpt = 1:length(tpt_array)
        iTpt
        ModelTript.parameters(9,:) = {'tpt',tpt_array(iTpt)};
    
        % Get FSP fit for bounds.
        ModelTript.solutionScheme = 'FSP';
        ModelTript.fspOptions.fspTol = 1e-6;
        ModelTript.fspOptions.bounds=[];
        [fspSolntpt,ModelTript.fspOptions.bounds] = ModelTript.solve;
    
        ModelTript.fspOptions.fspTol = inf;
        ModelTript.solutionScheme = 'fspSens';
        sensSoln_Tript{iTpt} = ModelTript.solve(fspSoln.stateSpace);
    
        % Solve for Fisher Information Matrix at all Time Points
        ModelTript.pdoOptions.unobservedSpecies = {'x1','x2'};
        fims_Tript = ModelTript.computeFIM(sensSoln_Tript{iTpt}.sens);
        FIM_Tript = ModelTript.evaluateExperiment(fims_Tript,ModelTript.dataSet.nCells);
        expectedDetCov(iTpt) = det(FIM_Tript^(-1));
    end
    
    % figure slide 27
    if showFigures
        figure;plot(tpt_array,expectedDetCov); hold on
        set(gca,'yscale','log')
        ylim = get(gca,'ylim');
        for it = 1:length(ModelTript.dataSet.times)
            plot(ModelTript.dataSet.times(it)*[1,1],ylim,'k--')
        end
    end
end

%% Save Data
if saveMHOutput % Save mhResults
    i=1;
    if exist(['EGRNT_mhast.mat'])
        while exist(['EGRNT_mhast','_v',num2str(i),'.mat'])
            i=i+1;
        end
    save(['EGRNT_mhast','_v',num2str(i),'.mat'],'mhResults');
    else
        save('EGRNT_mhast.mat','mhResults');
    end 
end

if saveMHOutput % Save Model
    j=1;
    if exist(['EGRNT_model.mat'])
        while exist(['EGRNT_model','_v',num2str(j),'.mat'])
            j=j+1;
        end
    save(['EGRNT_model','_v',num2str(j),'.mat'],'Model');
    else
        save('EGRNT_model','Model');
    end 
end

if saveSensOutput   % Save Sensitivity Results
    k=1;
    if exist(['EGRNT_sens.mat'])
        while exist(['EGRNT_sens','_v',num2str(k),'.mat'])
            k=k+1;
        end
    save(['EGRNT_sens','_v',num2str(k),'.mat'],'sensSoln');
    else
        save('EGRNT_sens.mat','sensSoln');
    end 
end

if saveSensOutput   % Save MLE Results
    n=1;
    if exist(['EGRNT_MLE.mat'])
        while exist(['EGRNT_MLE','_v',num2str(n),'.mat'])
            n=n+1;
        end
    save(['EGRNT_MLE','_v',num2str(n),'.mat'],'EGRNT_MLE');
    else
        save('EGRNT_MLE.mat','EGRNT_MLE');
    end 

end

if saveTriptOutput   % Save potimal Tript experiment Results
    n=1;
    if exist(['EGRNT_Tript_FIM.mat'])
        while exist(['EGRNT_Tript_FIM','_v',num2str(n),'.mat'])
            n=n+1;
        end
    save(['EGRNT_Tript_FIM','_v',num2str(n),'.mat'],'FIM_Tript');
    else
        save('EGRNT_Tript_FIM.mat','FIM_Tript');
    end 

end
