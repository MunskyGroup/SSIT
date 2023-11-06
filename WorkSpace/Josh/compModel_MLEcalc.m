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
    %% 
  starttime=-15000
  comp_Modelssa=comp_Model;
  comp_Modelssa.initialTime = starttime;
  comp_Modelssa.tSpan = [starttime,comp_Model.tSpan]
  comp_Modelssa.ssaOptions.useTimeVar=true
  comp_Modelssa.ssaOptions.signalUpdateRate=1
  comp_Modelssa.solutionScheme = 'SSA';  % Set solution scheme to SSA.
  comp_Modelssa.ssaOptions.Nexp = 200; 
  comp_Modelssa.ssaOptions.nSimsPerExpt = 50;
  % comp_Model.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
  comp_Modelssa.solve([],'CompDuspSSAData200Expts_B.csv'); 
  comp_Modelssa.solutionScheme = 'FSP';  % Set solution scheme back to FSP.

  %% Find MLE for each simulated data set.

  comp_Model.fittingOptions.modelVarsToFit = [1:4];
  nFrePars = length(comp_Model.fittingOptions.modelVarsToFit)
  MLE = zeros(1,nFrePars,comp_Modelssa.ssaOptions.Nexp);
  fMLE = inf(1,nFrePars,comp_Modelssa.ssaOptions.Nexp);
  %B = comp_Model.loadData('CompDuspSSAData200Expts.csv',{'x1','exp1_s1';'x2','exp1_s2';'x3','exp1_s3'});

for iExp = 1:comp_Modelssa.ssaOptions.Nexp
    iExp
    %GR&Gene&Dusp1   B = comp_Model.loadData('CompDuspSSAData200Expts.csv',{'x1',['exp',num2str(iExp),'_s1'];'x2',['exp',num2str(iExp),'_s2'];'x3',['exp',num2str(iExp),'_s3']}); % Link non-distorted data.
    %GR&Dusp1   B = comp_Model.loadData('CompDuspSSAData200Expts.csv',{'x1',['exp',num2str(iExp),'_s1'];'x3',['exp',num2str(iExp),'_s3']}); % Link non-distorted data.
    fitOptions = optimset;
    fitOptions.MaxIter = 100;
    fitOptions.Display = 'iter';

    B = comp_Model.loadData('CompDuspSSAData200Expts_B.csv',{'x3',['exp',num2str(iExp),'_s3']}); % Link non-distorted data.
    B.fittingOptions.timesToFit = [false,ones(1,length(comp_Model.tSpan),'logical')];
    B.tSpan = B.tSpan(2:end);
%     B.fittingOptions.modelVarsToFit = [2,3];
     if iExp==1
         x0 = [B.parameters{B.fittingOptions.modelVarsToFit,2}]';
     else
         x0 = squeeze(MLE(1,:,iExp-1)); 
     end
    [MLE(1,:,iExp),fMLE(1,:,iExp)] = B.maximizeLikelihood(x0,fitOptions);
end

%% Solve for Fisher Information Matrix at all Time Points
load MLE_result
B.fittingOptions.modelVarsToFit = [1:4];

B.solutionScheme = 'FSP';
B.fspOptions.fspTol = 1e-6;
B.fspOptions.bounds=[];
[fspSoln,B.fspOptions.bounds] = B.solve;

B.fspOptions.fspTol = inf;
B.solutionScheme = 'fspSens';
sensSoln = B.solve(fspSoln.stateSpace);

B.pdoOptions.unobservedSpecies = {'x1','x2'};
fims = B.computeFIM(sensSoln.sens);
FIM = B.evaluateExperiment(fims,B.dataSet.nCells);


%% Make Plots
fimFreePars = FIM(B.fittingOptions.modelVarsToFit,B.fittingOptions.modelVarsToFit);
B.makeMleFimPlot(squeeze(MLE(1,:,:)),fimFreePars,[4,3],0.95)