%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
%% Define SSIT Model
  Model = SSIT;
  Model.species = {'x1';'x2';'x3'};
  Model.initialCondition = [0;0;0];
  Model.propensityFunctions = {'kelt2*x2/(alpha+x2)*beta/(beta+x1)';'gelt2*x1';... % eLT-2
      'kelt7';'gelt7*x2';...% ELT-7
      'kgfp*x2/(alpha+x2)*beta/(beta+x1)';'ggfp*x3'}; %GFP
  Model.stoichiometry = [1,-1,0,0,0,0;...
                         0,0,1,-1,0,0;...
                         0,0,0,0,1,-1];
%   Model.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
  Model.parameters = ({'kelt2',1;...
      'kelt7',1;...
      'kgfp',1;...
      'gelt2',1;...
      'gelt7',1;...
      'ggfp',16;...
      'alpha',5; ...
      'beta',10});
  Model.fspOptions.initApproxSS = true;
  Model.tSpan=[0,100];

%% Solve the model using the FSP
  Model.solutionScheme = 'FSP';
  Model.fspOptions.fspTol = 1e-4;
  [fspSoln,Model.fspOptions.bounds] = Model.solve;

%% Load and Fit smFISH Data
  Model = Model.loadData('datasetDavidKingJune2023',{'x3','IntensityIntegers'},...
      {'genotype','ELT-2-GFP';'treatment','L4440'});
  Model.fspOptions.fspTol = inf;
  Model.fittingOptions.modelVarsToFit = 1:8;
  fitOptions = optimset('Display','iter','MaxIter',40);
  %%
  Model.fspOptions.bounds(6) = 40;

  Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(Model.maximizeLikelihood([],fitOptions));
  Model.makeFitPlot;



  %%
  ModelRNAi = Model;
  ModelRNAi.propensityFunctions = {'dRNAi*kelt2*x2/(alpha+x2)*beta/(beta+x1)';'gelt2*x1';... % eLT-2
      'kelt7';'gelt7*x2';...% ELT-7
      'kgfp*x2/(alpha+x2)*beta/(beta+x1)';'ggfp*x3'}; %GFP
  ModelRNAi.parameters(9,:) = {'dRNAi',0.1};
  ModelRNAi = ModelRNAi.loadData('datasetDavidKingJune2023',{'x3','IntensityIntegers'},...
      {'genotype','ELT-2-GFP';'treatment','ELT-2'});
  ModelRNAi.fittingOptions.modelVarsToFit = 7:9;
  ModelRNAi.fspOptions.bounds(4:6) = 30;
  ModelRNAi.parameters(ModelRNAi.fittingOptions.modelVarsToFit,2) = num2cell(ModelRNAi.maximizeLikelihood([],fitOptions));
  ModelRNAi.makeFitPlot;
  


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
