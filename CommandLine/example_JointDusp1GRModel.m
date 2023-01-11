%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
%% Define SSIT Model
  Model = SSIT;
  Model.species = {'x1';'x2';'x3'};  % GRnuc, geneOn, dusp1
  Model.initialCondition = [0;0;0];
  Model.propensityFunctions = {'(kcn0+kcn1*IDex)*(20-x1)';'knc*x1';...
      'kon*x1*(2-x2)';'koff*x2';'kr*x2';'gr*x3'};
  Model.stoichiometry = [ 1,-1, 0, 0, 0, 0;...
                          0, 0, 1,-1, 0, 0;...
                          0, 0, 0, 0, 1,-1];
  Model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'};
  Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
                       'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
  Model.fspOptions.initApproxSS = true;

  %% Solve the model using the FSP
  for i=1:10
      Model.solutionScheme = 'FSP';
      Model.fspOptions.fspTol = 1e-4;
      Model.fspOptions.verbose = 0;
      Model.fspOptions.bounds=[];
      [fspSoln,Model.fspOptions.bounds] = Model.solve;

      % Load and Fit smFISH Data
      %   Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'});
      Model.fspOptions.fspTol = inf;
      Model.fittingOptions.modelVarsToFit = 1:8;
      fitOptions = optimset('Display','iter','MaxIter',400);
      Model.parameters(1:8,2) = num2cell(Model.maximizeLikelihood([],fitOptions));
      Model.makeFitPlot;
  end

%% Metropolis Hastings to Quantify Parameter Uncertainty
  Model.fittingOptions.modelVarsToFit = 1:4;
  MHOptions = struct('numberOfSamples',1000,'burnin',500,'thin',3,...
                   'useFIMforMetHast',true,'suppressFSPExpansion',true);
  [~,~,mhResults] = Model.maximizeLikelihood([Model.parameters{1:4,2}]',...
                   MHOptions,'MetropolisHastings');
  Model.plotMHResults(mhResults);

%% Calibrate PDO from Multi-Modal Experimental Data
  ModelPDO = Model.calibratePDO('pdoCalibrationData.csv',...
    {'x2'},{'nTotal'},{'nType1'},'AffinePoiss',true);

%% Calculate CME Sensitivity to Parameter Variations
  Model.sensOptions.solutionMethod = 'finiteDifference';
  Model.solutionScheme = 'fspSens';
  Model.fspOptions.fspTol = 1e-6;
  [sensSoln,Model.fspOptions.bounds] = Model.solve;

%% Solve for Fisher Information Matrix at all Time Points
  Model.pdoOptions.unobservedSpecies = {'x1'};
  fims = Model.computeFIM(sensSoln.sens);
  FIM = Model.evaluateExperiment(fims,Model.dataSet.nCells);
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
 
