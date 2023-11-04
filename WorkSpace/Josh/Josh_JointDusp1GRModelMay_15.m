%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
%% Complex Dusp1 model
  Model = SSIT;
  Model.species = {'x1';'x2';'x3'};  % GRnuc, geneOn, dusp1
  Model.initialCondition = [0;0;0];
  Model.propensityFunctions = {'(kcn0+kcn1*IDex)';'knc*x1';...
      'kon*x1*(2-x2)';'koff*x2';'kr*x2';'gr*x3'};
  Model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'};
  Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
      'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
  Model.stoichiometry = [ 1,-1, 0, 0, 0, 0;...
                          0, 0, 1,-1, 0, 0;...
                          0, 0, 0, 0, 1,-1];
  Model.fspOptions.initApproxSS = true;

  %% Load saved mHast and Sensitivity
   mhResults = load('complex_dusp1_mhast.mat').mhResults;
   sensSoln = load('complex_dusp1_sens.mat').sensSoln;
   comp_Model = load('complex_dusp1_model.mat').Model;
%   %fspSoln = load('complex_dusp1_FSP.mat').fspSoln;
%   comp_Model.plotMHResults(mhResults);

  %mhResults = load('complex_dusp1_Trp_mhast.mat').mhResults;
  %sensSoln = load('complex_dusp1_Trp_sens.mat').sensSoln;
  %ModelTrypt = load('complex_dusp1_Trp_model.mat').ModelTrypt;
  
  %% Simutaneous Dusp1 + GR Measured
  comp_Model.pdoOptions.unobservedSpecies = {'x2'};
  fims_both = comp_Model.computeFIM(sensSoln.sens);
  FIM_both = comp_Model.evaluateExperiment(fims_both,comp_Model.dataSet.nCells);
  %comp_Model.plotMHResults(mhResults,FIM_both);
  %sgtitle('Simutaneous Dusp1 + GR Measured')

  nTotal_both = sum(comp_Model.dataSet.nCells);
  nCellsOpt_both = comp_Model.optimizeCellCounts(fims_both,nTotal_both,'TR[1:4]');% Case 1 TR[1:4], Case 2 [1:4]
  nCellsOpt_both_case2 = comp_Model.optimizeCellCounts(fims_both,nTotal_both,'[1:4]');% Case 1 TR[1:4], Case 2 [1:4]

  fimOpt_both = comp_Model.evaluateExperiment(fims_both,nCellsOpt_both);
  fimOpt_both_case2 = comp_Model.evaluateExperiment(fims_both,nCellsOpt_both_case2);

  comp_Model.plotMHResults(mhResults,{FIM_both,fimOpt_both});

    %% GR Measured 
  comp_Model.pdoOptions.unobservedSpecies = {'x3','x2'};
  fims_GR = comp_Model.computeFIM(sensSoln.sens);
  FIM_GR = comp_Model.evaluateExperiment(fims_GR,comp_Model.dataSet.nCells);
  %comp_Model.plotMHResults(mhResults,FIM_GR);
  %sgtitle('Dusp1 Measured')

  nTotal_GR = sum(comp_Model.dataSet.nCells);
  nCellsOpt_GR = comp_Model.optimizeCellCounts(fims_GR,nTotal_GR,'TR[1:4]');
  nCellsOpt_GR_case2 = comp_Model.optimizeCellCounts(fims_GR,nTotal_GR,'[1:4]');

  fimOpt_GR = comp_Model.evaluateExperiment(fims_GR,nCellsOpt_GR);
  fimOpt_GR_case2 = comp_Model.evaluateExperiment(fims_GR,nCellsOpt_GR_case2);

  comp_Model.plotMHResults(mhResults,{FIM_GR,fimOpt_GR});

  %% Dusp1 Measured
  comp_Model.pdoOptions.unobservedSpecies = {'x1','x2'};
  fims_dusp = comp_Model.computeFIM(sensSoln.sens);
  FIM_dusp = comp_Model.evaluateExperiment(fims_dusp,comp_Model.dataSet.nCells);
  %comp_Model.plotMHResults(mhResults,FIM_dusp);
  %gtitle('Dusp1 Measured')

  nTotal_dusp = sum(comp_Model.dataSet.nCells);
  nCellsOpt_dusp = comp_Model.optimizeCellCounts(fims_dusp,nTotal_dusp,'TR[1:4]');
  nCellsOpt_dusp_case2 = comp_Model.optimizeCellCounts(fims_dusp,nTotal_dusp,'[1:4]');

  fimOpt_dusp = comp_Model.evaluateExperiment(fims_dusp,nCellsOpt_dusp);
  fimOpt_dusp_case2 = comp_Model.evaluateExperiment(fims_dusp,nCellsOpt_dusp_case2);

  comp_Model.plotMHResults(mhResults,{FIM_dusp,fimOpt_dusp});

  %% Seperate dusp1 and GR measured (2X Cells)
  %   fims = comp_Model.computeFIM(sensSoln.sens);
  fims_dusp_and_gr_separate = [fims_dusp;fims_GR];
  
  % This is allowing us to use twice as many cells:
  FIM = FIM_GR +FIM_dusp;
  %comp_Model.plotMHResults(mhResults,FIM);
  %sgtitle('Seperate dusp1 and GR measured')

  % Here we will also allow twice as many cells to be measured with this
  % experiment in order to match the case above.
  nTotal = sum(comp_Model.dataSet.nCells);
  nCellsOpt_sep = comp_Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal*2,'TR[1:4]');
  nCellsOpt_sep_case2 = comp_Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal*2,'[1:4]');

  fimOpt = comp_Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep);
  fimOpt_case2 = comp_Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep_case2);

  comp_Model.plotMHResults(mhResults,{FIM,fimOpt});
 
  %% Seperate dusp1 and GR measured (1X Cells)
  %   fims = comp_Model.computeFIM(sensSoln.sens);
  fims_dusp_and_gr_separate = [fims_dusp;fims_GR];
  
  % This will fix the number of cells at half the original number for each
  % experiment.
  FIM = (FIM_GR+FIM_dusp)/2;
  %comp_Model.plotMHResults(mhResults,FIM);
  %sgtitle('Seperate dusp1 and GR measured')

  % Here we will also allow twice as many cells to be measured with this
  % experiment in order to match the case above.
  nTotal = sum(comp_Model.dataSet.nCells);
  nCellsOpt_sep = comp_Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal,'TR[1:4]');
  nCellsOpt_sep_case2 = comp_Model.optimizeCellCounts(fims_dusp_and_gr_separate,nTotal,'[1:4]');

  fimOpt = comp_Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep);
  fimOpt_case2 = comp_Model.evaluateExperiment(fims_dusp_and_gr_separate,nCellsOpt_sep_case2);

  comp_Model.plotMHResults(mhResults,{FIM,fimOpt});

  %% Dusp1 + GR + Gene state Measured
  comp_Model.pdoOptions.unobservedSpecies = {};
  fims_all = comp_Model.computeFIM(sensSoln.sens);
  FIM_all = comp_Model.evaluateExperiment(fims_all,comp_Model.dataSet.nCells);
  %comp_Model.plotMHResults(mhResults,FIM_all);
  %sgtitle('All Measured')

  nTotal_all = sum(comp_Model.dataSet.nCells);
  nCellsOpt_all = comp_Model.optimizeCellCounts(fims_all,nTotal,'TR[1:4]');
  nCellsOpt_all_case2 = comp_Model.optimizeCellCounts(fims_all,nTotal,'[1:4]');

  fimOpt_all = comp_Model.evaluateExperiment(fims_all,nCellsOpt_all);
  fimOpt_all_case2 = comp_Model.evaluateExperiment(fims_all,nCellsOpt_all_case2);

  comp_Model.plotMHResults(mhResults,{FIM_all,fimOpt_all});

  %%  Plot cell's per time stack
  timepoints = {'0min','10min','20min','30min','40min','50min','60min','75min','90min','120min','150min','180min'};

  figure()%case1
  designs = [comp_Model.dataSet.nCells;nCellsOpt_dusp;nCellsOpt_sep(1:12);nCellsOpt_both];
  designs_b = [comp_Model.dataSet.nCells;nCellsOpt_dusp;nCellsOpt_sep(1:12)+nCellsOpt_sep(13:24);nCellsOpt_both];


  g = bar(designs_b','barWidth',1);
  hold on
  h = bar(designs','barWidth',1);

  set(gca,'XTickLabel',timepoints,'FontSize',20)
  xlabel('Measurement timepoints')
  ylabel('Number of Cell')
  title('Case 1 (GR parameters known)')
  legend({'','',append('Seperate GR Experiments Optimized (',int2str(sum(nCellsOpt_sep(13:24))/sum(nCellsOpt_sep)*100),'% of cells)')...
      ,'','Equal cells at each timepoint','DUSP1 Optimized',append('Seperate DUSP1 Experiments Optimized (',int2str(sum(nCellsOpt_sep(1:12))/sum(nCellsOpt_sep)*100),'% of cells)'),'Simutaneous DUSP1/GR Experiments Optimized',''},'FontSize',20)

  figure()%case2
  designs = [comp_Model.dataSet.nCells;nCellsOpt_dusp_case2;nCellsOpt_sep_case2(1:12);nCellsOpt_both_case2];
  designs_b = [comp_Model.dataSet.nCells;nCellsOpt_dusp_case2;nCellsOpt_sep_case2(1:12)+nCellsOpt_sep_case2(13:24);nCellsOpt_both_case2];

  g = bar(designs_b','barWidth',1);
  hold on
  h = bar(designs','barWidth',1);

  set(gca,'XTickLabel',timepoints,'FontSize',20)
  xlabel('Measurement timepoints')
  ylabel('Number of Cell')
  title('Case 2 (GR parameters not known)')
  legend({'','',append('Seperate GR Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(13:24))/sum(nCellsOpt_sep_case2)*100),'% of cells)')...
      ,'','Equal cells at each timepoint','DUSP1 Optimized',append('Seperate DUSP1 Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(1:12))/sum(nCellsOpt_sep_case2)*100),'% of cells)'),'Simutaneous DUSP1/GR Experiments Optimized',''},'FontSize',20)

%% bar plots case 1
  y = [det(inv(FIM_dusp(1:4,1:4))),det(inv(fimOpt_dusp(1:4,1:4)));...
      det(FIM(1:4,1:4)),det(fimOpt(1:4,1:4)); det(FIM_both(1:4,1:4)), det(fimOpt_both(1:4,1:4))];
  x = categorical({'Dusp1', 'Seperate GR and Dusp1 experiment', 'Simultaneous GR and Dusp1 experiment'});
  bar(x,y,'barWidth',1)
  set(gca, 'YScale', 'log')%,'ylim',[1e-36,1e-20])
  ax=gca;
  ax.FontSize = 20;
  legend({'FIM','FIM optimized'},'FontSize', 24)
  ylabel('|FIM|')
  title('Case 1 (GR parameters know)')

  %% bar plots case 2
  ev = [eye(4),zeros(4)];
  y = [det(ev*inv(FIM_dusp(1:8,1:8))*ev'),det(ev*inv(fimOpt_dusp_case2(1:8,1:8))*ev');...
      det(ev*inv(FIM(1:8,1:8))*ev'),det(ev*inv(fimOpt_case2(1:8,1:8))*ev'); det(ev*inv(FIM_both(1:8,1:8))*ev'), det(ev*inv(fimOpt_both_case2(1:8,1:8))*ev')];
  x = categorical({'Dusp1', 'Seperate GR and Dusp1 experiment', 'Simultaneous GR and Dusp1 experiment'});
  bar(x,y,'barWidth',1)
  set(gca, 'YScale', 'log','ylim',[1e-36,1e-10])
  ax2=gca;
  ax2.FontSize = 25;
  legend({'FIM','FIM optimized'},'FontSize', 24)
  ylabel('|COV|')
  title('Case 2 (GR parameters not know')

  %% det(FIM^-1) vs # of cells 

  n = logspace(2,7,1000);
  for i=1:length(n)
     inv_FIM_case1(i) = det(inv((n(i)/nTotal_all)*FIM_dusp(1:4,1:4)));
     inv_FIM_case2(i) = det(ev*inv((n(i)/nTotal_all)*FIM_dusp(1:8,1:8))*ev');

  end
  lnwdt=2;
  figure()%Case 1
  plot(n,inv_FIM_case1,'k','LineWidth',lnwdt)
  yline(det(inv(FIM_dusp(1:4,1:4))),'b','LineWidth',lnwdt);
  yline(det(inv(FIM(1:4,1:4))),'g','LineWidth',lnwdt);
  yline(det(inv(FIM_both(1:4,1:4))),'m','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-35,1e-10],'xlim',[0,1e7],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'','DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment'});
  title('Case 1')

  figure()%Case 2
  plot(n,inv_FIM_case2,'k','LineWidth',lnwdt)
  yline(det(ev*inv(FIM_dusp(1:8,1:8))*ev'),'b','LineWidth',lnwdt);
  yline(det(ev*inv(FIM(1:8,1:8))*ev'),'g','LineWidth',lnwdt);
  yline(det(ev*inv(FIM_both(1:8,1:8))*ev'),'m','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-35,1e-10],'xlim',[0,1e7],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'','DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment'});
  title('Case 2')

  figure()%Case 1 opt
  plot(n,inv_FIM_case1,'k','LineWidth',lnwdt)
  yline(det(inv(FIM_dusp(1:4,1:4))),'b','LineWidth',lnwdt);
  yline(det(inv(FIM(1:4,1:4))),'g','LineWidth',lnwdt);
  yline(det(inv(FIM_both(1:4,1:4))),'m','LineWidth',lnwdt);
  yline(det(inv(fimOpt_dusp(1:4,1:4))),'--b','LineWidth',lnwdt);
  yline(det(inv(fimOpt(1:4,1:4))),'--g','LineWidth',lnwdt);
  yline(det(inv(fimOpt_both(1:4,1:4))),'--m','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-35,1e-10],'xlim',[0,1e7],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'','DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment'});
  title('Case 1')

  figure()%Case 2 opt
  plot(n,inv_FIM_case2,'k','LineWidth',lnwdt)
  yline(det(ev*inv(FIM_dusp(1:8,1:8))*ev'),'b','LineWidth',lnwdt);
  yline(det(ev*inv(FIM(1:8,1:8))*ev'),'g','LineWidth',lnwdt);
  yline(det(ev*inv(FIM_both(1:8,1:8))*ev'),'m','LineWidth',lnwdt);
  yline(det(ev*inv(fimOpt_dusp_case2(1:8,1:8))*ev'),'--b','LineWidth',lnwdt);
  yline(det(ev*inv(fimOpt_case2(1:8,1:8))*ev'),'--g','LineWidth',lnwdt);
  yline(det(ev*inv(fimOpt_both_case2(1:8,1:8))*ev'),'--m','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-35,1e-10],'xlim',[0,1e7],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'','DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment'});
  title('Case 2')

%% det(FIM^-1) vs # of cells no stack
 
% Plot the inv_FIM_case1 of each experiment with horizontial line at 1e-20
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
  figure()%Case 1
  plot(n,inv_FIM_case1_dusp,'b',n,inv_FIM_case1_sep,'g',n,inv_FIM_case1_both,'m','LineWidth',lnwdt)
  yline(9.97e-25,'k','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-35,1e-10],'xlim',[0,1e8],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment','Desired Information for Experiments'},'FontSize',20);
  title('Case 1')

  figure()%Case 2
  plot(n,inv_FIM_case2_dusp,'b',n,inv_FIM_case2_sep,'g',n,inv_FIM_case2_both,'m','LineWidth',lnwdt)
  yline(1e-14,'k','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-20,1e-10],'xlim',[0,1e8],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment','Desired Information for Experiments'},'FontSize',20);
  title('Case 2')

  figure()%Case 1 opt
  plot(n,inv_FIM_case1_dusp,'b',n,inv_FIM_case1_sep,'g',n,inv_FIM_case1_both,'m',...
      n,inv_FIM_case1_dusp_opt,'--b',n,inv_FIM_case1_sep_opt,'--g',n,inv_FIM_case1_both_opt,'--m','LineWidth',lnwdt)
  yline(9.97e-25,'k','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-35,1e-10],'xlim',[0,1e8],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment','','','','Desired Information for Experiments'},'FontSize',20);
  title('Case 1')

  figure()%Case 2 opt
  plot(n,inv_FIM_case2_dusp,'b',n,inv_FIM_case2_sep,'g',n,inv_FIM_case2_both,'m',...
      n,inv_FIM_case2_dusp_opt,'--b',n,inv_FIM_case2_sep_opt,'--g',n,inv_FIM_case2_both_opt,'--m','LineWidth',lnwdt)
  yline(1e-14,'k','LineWidth',lnwdt);
  set(gca, 'YScale', 'log','XScale', 'log','ylim',[1e-20,1e-10],'xlim',[0,1e8],'FontSize', 24)
  xlabel('Number of Cells Measured')
  ylabel('|COV|')
  legend({'DUSP1 experiment','Seperate GR and Dusp1 experiment','Simultaneous GR and Dusp1 experiment','','','','Desired Information for Experiments'},'FontSize',20);
  title('Case 2')
  
  %% |COV| plotted vs time
  trpt_times_vect = linspace(0,180,10);
  for i = length(trpt_times_vect):-1:1
      i
      tic
      ModelTrypt_t{i} = ModelTrypt;
      ModelTrypt_t{i}.parameters{9,2} = trpt_times_vect(i);
      ModelTrypt_t{i}.sensOptions.solutionMethod = 'finiteDifference';
      ModelTrypt_t{i}.solutionScheme = 'fspSens';
      ModelTrypt_t{i}.fspOptions.fspTol = 1e-6;
      [sensSoln_t{i},ModelTrypt_t{i}.fspOptions.bounds] = ModelTrypt_t{i}.solve;
      %fimResults{i} = ModelTrypt_t{i}.computeFIM(sensSoln{i});
      fimResults{i} = ModelTrypt_t{i}.computeFIM(sensSoln.sens);
      FIM_all{i} = ModelTrypt_t{i}.evaluateExperiment(fimResults{i},ModelTrypt_t{i}.dataSet.nCells);
      FIM_det(i) = det(FIM_all{i}(1:8,1:8))
      
      toc
  end
  plot(trpt_times_vect,FIM_det)
  %%
  set(gca,'XTickLabel',timepoints,'FontSize',20)
  xlabel('Measurement timepoints')
  ylabel('Number of Cell')
  title('Case 2 (GR parameters not known)')
  legend({'','',append('Seperate GR Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(13:24))/sum(nCellsOpt_sep_case2)*100),'% of cells)')...
      ,'','Equal cells at each timepoint','DUSP1 Optimized',append('Seperate DUSP1 Experiments Optimized (',int2str(sum(nCellsOpt_sep_case2(1:12))/sum(nCellsOpt_sep_case2)*100),'% of cells)'),'Simutaneous DUSP1/GR Experiments Optimized',''},'FontSize',20)
  %% 8 parameter mhast
  mhResults = load('8_par_complex_dusp1_mhast.mat').mhResults;
  sensSoln = load('8_par_complex_dusp1_sens.mat').sensSoln;
  par_8_Model = load('8_par_complex_dusp1_model.mat').Model;

  par_8_Model.plotMHResults(mhResults);
  
  %% Solve the model using the FSP
  Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'});
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
      fitOptions = optimset('Display','iter','MaxIter',10);
      Model.parameters(1:8,2) = num2cell(Model.maximizeLikelihood([],fitOptions));
      Model.makeFitPlot;
  end

%% Metropolis Hastings to Quantify Parameter Uncertainty
  Model.fittingOptions.modelVarsToFit = 1:8;
  for i =1:5
      MHOptions = struct('numberOfSamples',200,'burnin',0,'thin',1,...
          'useFIMforMetHast',true,'suppressFSPExpansion',true);
      [bestParsFound,~,mhResults] = Model.maximizeLikelihood([Model.parameters{Model.fittingOptions.modelVarsToFit,2}]',...
          MHOptions,'MetropolisHastings');
      Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
      Model.plotMHResults(mhResults);
  end
  

%% Calibrate PDO from Multi-Modal Experimental Data
  ModelPDO = Model.calibratePDO('pdoCalibrationData.csv',...
    {'x2'},{'nTotal'},{'nType1'},'AffinePoiss',true);

%% Calculate CME Sensitivity to Parameter Variations
  Model.sensOptions.solutionMethod = 'finiteDifference';
  Model.solutionScheme = 'fspSens';
  Model.fspOptions.fspTol = 1e-6;
  [sensSoln,Model.fspOptions.bounds] = Model.solve;

%% Solve for Fisher Information Matrix at all Time Points
  Model.pdoOptions.unobservedSpecies = {'x1','x2'};
  fims = Model.computeFIM(sensSoln.sens);
%   numberCellsPerTimePoint = Model.dataSet.nCells;
  numberCellsPerTimePoint 
  FIM = Model.evaluateExperiment(fims,numberCellsPerTimePoint);
  Model.plotMHResults(mhResults,FIM);

%% Optimize Experiment Design (Same Number of Cells and Timepoints)
  nTotal = sum(Model.dataSet.nCells);
  nCellsOpt = Model.optimizeCellCounts(fims,nTotal,'TR[1:4]');
  fimOpt = Model.evaluateExperiment(fims,nCellsOpt);
  Model.plotMHResults(mhResults,{FIM,fimOpt});
%%
  close all
  bar([1:12],AAA(1,:),.45,'k'); hold on
  bar([1:12]+0.45,AAA(2,:),.45,'c'); hold on
  set(gca,'xtick',[1:12]+0.225,'XTickLabel',Model.tSpan,'fontsize',15)
  legend('Intuitive Design','Optimized Design')

%% Tryptolide Experiment

% Complex Dusp1 model
  Model = SSIT;
  Model.species = {'x1';'x2';'x3'};  % GRnuc, geneOn, dusp1
  Model.initialCondition = [0;0;0];
  Model.propensityFunctions = {'(kcn0+kcn1*IDex)';'knc*x1';...
      'kon*x1*(2-x2)';'koff*x2';'kr*x2';'gr*x3'};
  Model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'};
  Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
      'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
  Model.stoichiometry = [ 1,-1, 0, 0, 0, 0;...
      0, 0, 1,-1, 0, 0;...
      0, 0, 0, 0, 1,-1];
  Model.fspOptions.initApproxSS = true;

%%  Setting data
ModelTrypt = Model;
% load dataset
ModelTrypt = ModelTrypt.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'});
ModelTrypt.fittingOptions.modelVarsToFit = 1:9;

ModelTrypt.propensityFunctions(5) = {'kr*x2*Itrypt'};
ModelTrypt.inputExpressions(2,:) = {'Itrypt','(t<tpt)'};
tpt_array = 20:20:180;
ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';

%% Running FSP fits
for i=1:5
    ModelTrypt.solutionScheme = 'FSP';
    ModelTrypt.fspOptions.fspTol = 1e-6;
    ModelTrypt.parameters(9,:) = {'tpt',180};
    ModelTrypt.fspOptions.bounds=[];
    [fspSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;
    ModelTrypt.fspOptions.bounds

% Load and Fit smFISH Data
      ModelTrypt.fspOptions.fspTol = inf;
      ModelTrypt.fittingOptions.modelVarsToFit = 1:9;
      fitOptions = optimset('Display','iter','MaxIter',200);
      ModelTrypt.parameters(1:9,2) = num2cell(ModelTrypt.maximizeLikelihood([],fitOptions));
      ModelTrypt.makeFitPlot;
end


%% Metropolis Hastings to Quantify Parameter Uncertainty (TRP)
  ModelTrypt.fittingOptions.modelVarsToFit = 1:9;
  for i =1:5
      MHOptions = struct('numberOfSamples',1000,'burnin',0,'thin',1,...
          'useFIMforMetHast',true,'suppressFSPExpansion',true);
      [bestParsFound,~,mhResults] = ModelTrypt.maximizeLikelihood([ModelTrypt.parameters{ModelTrypt.fittingOptions.modelVarsToFit,2}]',...
          MHOptions,'MetropolisHastings');
      ModelTrypt.parameters(ModelTrypt.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
      ModelTrypt.plotMHResults(mhResults);
  end
  

%% Calibrate PDO from Multi-Modal Experimental Data (TRP)
  ModelTryptPDO = ModelTrypt.calibratePDO('pdoCalibrationData.csv',...
    {'x2'},{'nTotal'},{'nType1'},'AffinePoiss',true);

%% Calculate CME Sensitivity to Parameter Variations (TRP)
  ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';
  ModelTrypt.solutionScheme = 'fspSens';
  ModelTrypt.fspOptions.fspTol = 1e-6;
  [sensSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;










%% fit not needed yet?

ModelTrypt.solutionScheme = 'fspSens';
ModelTrypt.fspOptions.fspTol = 9e-5;
ModelTrypt.pdoOptions.unobservedSpecies = {'x1','x2'};
for itpt = 1:length(tpt_array)
  ModelTrypt.parameters(9,:) = {'tpt',tpt_array(itpt)};
  [sensSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve(fspSoln.stateSpace);
  fims = ModelTrypt.computeFIM(sensSoln.sens);
  FIM = ModelTrypt.evaluateExperiment(fims,ModelTrypt.dataSet.nCells);
  fimChosenPars = FIM(ModelTrypt.fittingOptions.modelVarsToFit,...
      ModelTrypt.fittingOptions.modelVarsToFit);
  exptValue(itpt) = det(fimChosenPars)
end

%% hmm...
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
      fitOptions = optimset('Display','iter','MaxIter',10);
      Model.parameters(1:8,2) = num2cell(Model.maximizeLikelihood([],fitOptions));
      Model.makeFitPlot;
  end
 
