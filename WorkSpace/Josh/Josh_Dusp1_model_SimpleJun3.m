%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
  addpath('../')
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


  simple_Model = Model;
  tpt_array = 20:20:180;
  simple_Model.tSpan = tpt_array;

  %% Load saved mHast and Sensitivity
  %mhResults = load('simple_dusp1_mhast.mat').mhResults;
  %sensSoln = load('simple_dusp1_sens.mat').sensSoln;
  %simple_Model = load('simple_dusp1_model.mat').Model;
  %fspSoln = load('complex_dusp1_FSP.mat').fspSoln;
  %simple_Model.plotMHResults(mhResults);

  %% Define Prior
  % Define Prior on Model Parameters (log10-normal)
%   Model.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01;...
%                        'a1',0.4;'r1',0.04;'r2',0.1});
  muLog10Prior = [-1,-1,1,-2,-1,-2,-1];
  sigLog10Prior = [2 2 1 1 2 2 2];
  simple_Model.fittingOptions.modelVarsToFit = 1:7;
  muLog10Prior = muLog10Prior(simple_Model.fittingOptions.modelVarsToFit);
  sigLog10Prior = sigLog10Prior(simple_Model.fittingOptions.modelVarsToFit);
  simple_Model.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
  
  %% Solve the model using the FSP
  simple_Model = simple_Model.loadData('../../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});  
  simple_Model.fittingOptions.modelVarsToFit = 1:7;

  for i=1:5
      simple_Model.solutionScheme = 'FSP';
      simple_Model.fspOptions.fspTol = 1e-4;
      simple_Model.fspOptions.verbose = 0;
      simple_Model.fspOptions.bounds=[];
      [fspSoln,simple_Model.fspOptions.bounds] = simple_Model.solve;
      simple_Model.fspOptions.bounds

      % Load and Fit smFISH Data
      simple_Model.fspOptions.fspTol = inf;
      fitOptions = optimset('Display','iter','MaxIter',100);
      simple_Model.parameters(simple_Model.fittingOptions.modelVarsToFit,2) =...
          num2cell(simple_Model.maximizeLikelihood(...
          [simple_Model.parameters{simple_Model.fittingOptions.modelVarsToFit,2}],...
          fitOptions));
  end
  simple_Model.makeFitPlot;

%% Metropolis Hastings to Quantify Parameter Uncertainty
  simple_Model.fittingOptions.modelVarsToFit = 1:7;
%   for i =1:5
      MHOptions = struct('numberOfSamples',20,'burnin',0,'thin',1,...
          'useFIMforMetHast',true,'suppressFSPExpansion',true);
      [bestParsFound,~,mhResults] = simple_Model.maximizeLikelihood([simple_Model.parameters{simple_Model.fittingOptions.modelVarsToFit,2}]',...
          MHOptions,'MetropolisHastings');
      simple_Model.parameters(simple_Model.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
      simple_Model.plotMHResults(mhResults);
%   end
  

%% Calibrate PDO from Multi-Modal Experimental Data
  ModelPDO = simple_Model.calibratePDO('pdoCalibrationData.csv',...
    {'x2'},{'nTotal'},{'nType1'},'AffinePoiss',true);

%% Calculate CME Sensitivity to Parameter Variations
  simple_Model.sensOptions.solutionMethod = 'finiteDifference';
  simple_Model.solutionScheme = 'fspSens';
  simple_Model.fspOptions.fspTol = 1e-6;
  [sensSoln,simple_Model.fspOptions.bounds] = simple_Model.solve;

%% Solve for Fisher Information Matrix at all Time Points
  simple_Model.pdoOptions.unobservedSpecies = {};
  fims = simple_Model.computeFIM(sensSoln.sens);
  FIM = simple_Model.evaluateExperiment(fims,simple_Model.dataSet.nCells);
  simple_Model.plotMHResults(mhResults,FIM);

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
ModelTrypt = simple_Model;
ModelTrypt.propensityFunctions(3) = {'kr*x1*Itrypt'};
ModelTrypt.inputExpressions(2,:) = {'Itrypt','(t<tpt)'};


ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';

tpt_array = 30:30:180;
sensSoln = cell(1,length(tpt_array));
for iTpt = 1:length(tpt_array)
    ModelTrypt.parameters(8,:) = {'tpt',tpt_array(iTpt)};

    % Get FSP fit for bounds.
    ModelTrypt.solutionScheme = 'FSP';
    ModelTrypt.fspOptions.fspTol = 1e-6;
    ModelTrypt.fspOptions.bounds=[];
    [fspSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;

    ModelTrypt.fspOptions.fspTol = inf;
    ModelTrypt.solutionScheme = 'fspSens';
    sensSoln{iTpt} = ModelTrypt.solve(fspSoln.stateSpace);

    %% Solve for Fisher Information Matrix at all Time Points
    ModelTrypt.pdoOptions.unobservedSpecies = {'x1'};
    fims = ModelTrypt.computeFIM(sensSoln{iTpt}.sens);
    FIM = ModelTrypt.evaluateExperiment(fims,ModelTrypt.dataSet.nCells);

    expectedDetCov(iTpt) = det(FIM^(-1))
end
figure;plot(tpt_array,expectedDetCov); hold on
set(gca,'yscale','log')
ylim = get(gca,'ylim');
for it = 1:length(ModelTrypt.dataSet.times)
    plot(ModelTrypt.dataSet.times(it)*[1,1],ylim,'k--')
end

%%  Setting data
ModelTrypt = Model;
% load dataset
ModelTrypt = ModelTrypt.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'});
ModelTrypt.fittingOptions.modelVarsToFit = 1:8;

ModelTrypt.propensityFunctions(3) = {'kr*x1*Itrypt'};
ModelTrypt.inputExpressions(2,:) = {'Itrypt','(t<tpt)'};
tpt_array = 20:20:180;
ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';

%% Running FSP fits
for i=1:5
    ModelTrypt.solutionScheme = 'FSP';
    ModelTrypt.fspOptions.fspTol = 1e-6;
    ModelTrypt.parameters(8,:) = {'tpt',180};
    ModelTrypt.fspOptions.bounds=[];
    [fspSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;
    ModelTrypt.fspOptions.bounds

% for itpt = 1:length(tpt_array)
%   ModelTrypt.parameters(8,:) = {'tpt',tpt_array(itpt)};
%   [sensSoln{itpt},ModelTrypt.fspOptions.bounds] = ModelTrypt.solve(fspSoln.stateSpace);
%   fims{itpt} = ModelTrypt.computeFIM(sensSoln{itpt}.sens);
%   FIM = ModelTrypt.evaluateExperiment(fims{itpt},ModelTrypt.dataSet.nCells);
%   fimChosenPars = FIM(ModelTrypt.fittingOptions.modelVarsToFit,...
%       ModelTrypt.fittingOptions.modelVarsToFit);
%   exptValue(itpt) = det(fimChosenPars)
% end


end


%% Metropolis Hastings to Quantify Parameter Uncertainty (TRP)
  ModelTrypt.fittingOptions.modelVarsToFit = 1:8;
  for i =1:5
      MHOptions = struct('numberOfSamples',1000,'burnin',0,'thin',1,...
          'useFIMforMetHast',true,'suppressFSPExpansion',true);
      [bestParsFound,~,mhResults] = ModelTrypt.maximizeLikelihood([ModelTrypt.parameters{ModelTrypt.fittingOptions.modelVarsToFit,2}]',...
          MHOptions,'MetropolisHastings');
      ModelTrypt.parameters(ModelTrypt.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
      ModelTrypt.plotMHResults(mhResults);
  end


%%
[sensSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;


%% save data
save('simple_dusp1_Trp_model.mat','ModelTrypt')
save('simple_dusp1_Trp_fspSoln.mat','fspSoln')
save('simple_dusp1_Trp_sens.mat','sensSoln')
save('simple_dusp1_Trp_mhast.mat','mhResults')

% %%
% ModelTrypt.solutionScheme = 'fspSens';
% ModelTrypt.fspOptions.fspTol = 9e-5;
% for itpt = 1:length(tpt_array)
%   ModelTrypt.parameters(8,:) = {'tpt',tpt_array(itpt)};
%   [sensSoln{itpt},ModelTrypt.fspOptions.bounds] = ModelTrypt.solve(fspSoln.stateSpace);
%   fims{itpt} = ModelTrypt.computeFIM(sensSoln{itpt}.sens);
%   FIM = ModelTrypt.evaluateExperiment(fims{itpt},ModelTrypt.dataSet.nCells);
%   fimChosenPars = FIM(ModelTrypt.fittingOptions.modelVarsToFit,...
%       ModelTrypt.fittingOptions.modelVarsToFit);
%   exptValue(itpt) = det(fimChosenPars)
% end
% 
%  %plot det*(FIM-1) vs tpt_array
