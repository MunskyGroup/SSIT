%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
%% Define SSIT Model
  Model = SSIT;
  Model.species = {'x1';'x2'};
  Model.initialCondition = [0;0];
  Model.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
  Model.stoichiometry = [1,-1,0,0;0,0,1,-1];
  Model.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
  Model.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01;...
                       'a1',0.4;'r1',0.04;'r2',0.1});
  Model.fspOptions.initApproxSS = true;

%% Solve the model using the FSP
  Model.solutionScheme = 'FSP';
  Model.fspOptions.fspTol = 1e-4;
  [fspSoln,Model.fspOptions.bounds] = Model.solve;

%% Load and Fit smFISH Data
  Model = Model.loadData('../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});
  Model.fspOptions.fspTol = inf;
  Model.fittingOptions.modelVarsToFit = 1:7;
  fitOptions = optimset('Display','iter','MaxIter',4000);
  %%
  Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(Model.maximizeLikelihood([],fitOptions));
  Model.makeFitPlot;

%% Metropolis Hastings to Quantify Parameter Uncertainty
  Model.fittingOptions.modelVarsToFit = 1:4;
  Model.sensOptions.solutionMethod = 'finiteDifference';
  MHOptions = struct('numberOfSamples',1000,'burnin',500,'thin',3,...
                   'useFIMforMetHast',true,'suppressFSPExpansion',true,'CovFIMscale',.6);
  [~,~,mhResults] = Model.maximizeLikelihood([Model.parameters{1:4,2}]',...
                   MHOptions,'MetropolisHastings');
  Model.plotMHResults(mhResults);

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
 
%% Calibrate PDO from Multi-Modal Experimental Data
  ModelPDOSpots = Model.calibratePDO('ExampleDataSets/pdoCalibrationData.csv',...
    {'x2'},{'nTotal'},{'nSpots0'},'AffinePoiss',true);

  ModelPDOIntens = Model.calibratePDO('ExampleDataSets/pdoCalibrationData.csv',...
    {'x2'},{'nTotal'},{'intens1'},'AffinePoiss',true,[1,1000,1]);

%%
fimsPDOSpot = ModelPDOSpots.computeFIM(sensSoln.sens);
fimPDOSpots = ModelPDOSpots.evaluateExperiment(fimsPDOSpot,nCellsOpt);
Model.plotMHResults(mhResults,{FIM,fimOpt,fimPDOSpots});

fimsPDOIntens = ModelPDOIntens.computeFIM(sensSoln.sens);
fimPDOIntens = ModelPDOIntens.evaluateExperiment(fimsPDOIntens,nCellsOpt);
Model.plotMHResults(mhResults,{FIM,fimOpt,fimPDOSpots,fimPDOIntens});

fimsPDOIntens = ModelPDOIntens.computeFIM(sensSoln.sens);
fimPDOIntens2x = ModelPDOIntens.evaluateExperiment(fimsPDOIntens,2.218*nCellsOpt);
Model.plotMHResults(mhResults,{FIM,fimOpt,fimPDOSpots,fimPDOIntens,fimPDOIntens2x});
