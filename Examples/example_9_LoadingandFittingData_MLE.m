%% SSIT/Examples/example_9_LoadingandFittingData_MLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.3: Loading and fitting time-varying STL1 yeast data 
%   * Maximize the likelihood L(D|theta) and use the maximum likelihood
%     estimate (MLE) to fit the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, and loaded data from 
% example_8_LoadingandFittingData_DataLoading
% clear
% close all
addpath(genpath('../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading

%% Load pre-computed FSP solutions & loaded data:
% load('example_8_LoadingandFittingData.mat')
 
% View model summary:
% Model_data.summarizeModel
% STL1_data.summarizeModel
STL1_4state_data.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit experimental data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make new copies of our model:
% Model_MLE = Model_data;
% STL1_MLE = STL1_data;
STL1_4state_MLE = STL1_4state_data;

% Specify how many model parameters will be fit (the rest will be fixed):
fitpars = 13;
STL1_4state_MLE.fittingOptions.modelVarsToFit = 1:fitpars;

% Set fitOptions, with the maximum allowable number of iterations to fit:
fitOptions = optimset('Display','iter','MaxIter',3000);

% Store parameters for fitting:
% Model_pars = cell2mat(Model_MLE.parameters(1:4,2));
% STL1_pars = cell2mat(STL1_MLE.parameters(1:8,2));
STL1_4state_pars = cell2mat(STL1_4state_MLE.parameters(1:fitpars,2));

%% Compute the MLEs:
% [Model_pars,Model_likelihood] = ...
%  Model_MLE.maximizeLikelihood(Model_pars,fitOptions);
% 
% [STL1_pars,STL1_likelihood] = ...
%  STL1_MLE.maximizeLikelihood(STL1_pars,fitOptions);

[STL1_4state_pars,STL1_4state_likelihood] = ...
 STL1_4state_MLE.maximizeLikelihood(STL1_4state_pars,fitOptions);
% Note: Should see an MLE of -21961.2 at the end

% Update parameters:
% for j=1:length(Model_pars)
%     Model_MLE.parameters{j,2} = Model_pars(j);
% end
% 
% for k=1:length(STL1_pars)
%     STL1_MLE.parameters{k,2} = STL1_pars(k);
% end

for l=1:length(STL1_4state_pars)
    STL1_4state_MLE.parameters{l,2} = STL1_4state_pars(l);
end

% Make plots of the parameter fits from the MLEs:
% Model_MLE.plotFits([], "all", [], {'linewidth',2},...
%     Title='Bursting Gene', YLabel='Molecule Count',...
%     LegendLocation='northeast', LegendFontSize=12);
% 
% STL1_MLE.plotFits([], "all", [], {'linewidth',2},...
%     Title='STL1', YLabel='Molecule Count',...
%     LegendLocation='northeast', LegendFontSize=12);

STL1_4state_MLE.plotFits([],"all",[],{'linewidth',2},TitleFontSize=24,...
    Title='4-state STL1 (MLE)', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=18, AxisLabelSize=20);

%% Save models & MLEs:
saveNames = unique({'Model_MLE'
    'STL1_MLE'
    'STL1_4state_MLE'
    'Model_pars'
    'STL1_pars'
    'STL1_4state_pars'
    'Model_likelihood'
    'STL1_likelihood'
    'STL1_4state_likelihood'
    });
    
save('example_9_LoadingandFittingData_MLE',saveNames{:})