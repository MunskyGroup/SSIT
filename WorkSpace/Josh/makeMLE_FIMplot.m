%% Clear work space
clc, clear

%% Load data for plotting MLE results
MLE = load('MLE_result.mat').MLE;
sensSoln = load('complex_dusp1_sens.mat').sensSoln;
comp_Model = load('complex_dusp1_model.mat').Model;

%% Making MLE FIM plot with MH points
fims = comp_Model.computeFIM(sensSoln.sens);
[FIM,sFIMcov,fimMetrics] = comp_Model.evaluateExperiment(fims,300);
fimFreePars = FIM(comp_Model.fittingOptions.modelVarsToFit,comp_Model.fittingOptions.modelVarsToFit);
comp_Model.makeMleFimPlot(squeeze(MLE(1,:,:)),fimFreePars,[4,3],0.95)