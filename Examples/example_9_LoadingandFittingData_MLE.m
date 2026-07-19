%% SSIT/Examples/example_9_LoadingandFittingData_MLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.3.3: Loading and fitting time-varying STL1 yeast data 
%   * Maximize the likelihood L(D|theta) and use the maximum likelihood
%     estimate (MLE) to fit the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, and loaded data from 
% example_8_LoadingandFittingData_DataLoading
% clear
% close all

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading

%% Load pre-computed FSP solutions & loaded data:
load('example_8_LoadingandFittingData_DataLoading.mat')
 
% View model summary:
% Model.summarizeModel
% STL1.summarizeModel
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit experimental data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify how many model parameters will be fit (the rest will be fixed):
fitpars = 13;
STL1_4state.fittingOptions.modelVarsToFit = 1:fitpars;

% Set fitOptions, with the maximum allowable number of iterations to fit:
fitOptions = optimset('Display','iter','MaxIter',2000);

%% Compute the MLEs:
% Model =Model.maximizeLikelihood(fitOptions=fitOptions);
% STL1 = STL1.maximizeLikelihood(fitOptions=fitOptions);

STL1_4state = STL1_4state.maximizeLikelihood(fitOptions=fitOptions);
% Note: Should see an MLE of -24863.7 at the end

% Make plots of the parameter fits from the MLEs:
% Model.plotFits(plotType="all", lineProps={'linewidth',2},...
%     Title='Bursting Gene', YLabel='Molecule Count',...
%     LegendLocation='northeast', LegendFontSize=12);
% 
% STL1.plotFits(plotType="all", lineProps={'linewidth',2},...
%     Title='STL1', YLabel='Molecule Count',...
%     LegendLocation='northeast', LegendFontSize=12);

STL1_4state.plotFits(plotType="all", lineProps={'linewidth',2},...
    TitleFontSize=24, Title='4-state STL1 (MLE)', LegendFontSize=18,...
    YLabel='Molecule Count', LegendLocation='northeast', AxisLabelSize=20);

%% Save models & MLEs:
saveNames = unique({
    'STL1_4state'
    });
    
save('example_9_LoadingandFittingData_MLE',saveNames{:})