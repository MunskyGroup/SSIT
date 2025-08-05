%% example_9_LoadingandFittingData_MLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%   * Maximize the likelihood L(D|theta) and use the maximum likelihood
%     estimate (MLE) to fit the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions from 
% example_4_SolveSSITModels_FSP, and loaded data from 
% example_8_LoadingandFittingData_DataLoading
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading

% View model summary:
Model_data.summarizeModel
STL1_data.summarizeModel
STL1_data_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit experimental data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make new copies of our model, a simple bursting gene model and a more 
% complex bursting gene model that has a time-varying input signal that 
% turns on the STL1 gene:
Model_MLE = Model_data;
STL1_MLE = STL1_data;
STL1_MLE_4state = STL1_data_4state;
% Let's see which model better fits our data...

% Set fitOptions, with the maximum allowable number of iterations to fit:
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
Modelpars = cell2mat(Model_MLE.parameters(1:4,2));
STL1pars = cell2mat(STL1_MLE.parameters(1:7,2));
STL1pars_4state = cell2mat(STL1_MLE_4state.parameters(1:15,2));

%% Compute the MLEs:
[Modelpars,Model_likelihood] = ...
 Model_MLE.maximizeLikelihood(Modelpars,fitOptions);

[STL1pars,STL1_likelihood] = ...
 STL1_MLE.maximizeLikelihood(STL1pars,fitOptions);

[STL1pars_4state,STL1_likelihood_4state] = ...
 STL1_MLE_4state.maximizeLikelihood(STL1pars_4state,fitOptions);

% Update parameters:
for j=1:length(Modelpars)
    Model_MLE.parameters{j,2} = Modelpars(j);
end

for k=1:length(STL1pars)
    STL1_MLE.parameters{k,2} = STL1pars(k);
end

for l=1:length(STL1pars_4state)
    STL1_MLE_4state.parameters{l,2} = STL1pars_4state(l);
end

% Make plots of the parameter fits from the MLEs:
Model_MLE.makeFitPlot
STL1_MLE.makeFitPlot
STL1_MLE_4state.makeFitPlot