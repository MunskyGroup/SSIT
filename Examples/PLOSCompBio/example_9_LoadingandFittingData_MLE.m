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
example_8_LoadingandFittingData_DataLoading

%% Load pre-computed FSP solutions & loaded data:
load('example_8_LoadingandFittingData.mat')

% View model summary:
Model_data.summarizeModel
STL1_data.summarizeModel
STL1_4state_data.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit experimental data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make new copies of our model, a simple bursting gene model and a more 
% complex bursting gene model that has a time-varying input signal that 
% turns on the STL1 gene:
Model_MLE = Model_data;
STL1_MLE = STL1_data;
STL1_4state_MLE = STL1_4state_data;
% Let's see which model better fits our data...

% Set fitOptions, with the maximum allowable number of iterations to fit:
fitOptions = optimset('Display','iter','MaxIter',3000);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double

Model_pars = cell2mat(Model_MLE.parameters(1:4,2));
STL1_pars = cell2mat(STL1_MLE.parameters(1:7,2));
STL1_4state_pars = cell2mat(STL1_4state_MLE.parameters(1:15,2));

%% Compute the MLEs:
[Model_pars,Model_likelihood] = ...
 Model_MLE.maximizeLikelihood(Model_pars,fitOptions);

[STL1_pars,STL1_likelihood] = ...
 STL1_MLE.maximizeLikelihood(STL1_pars,fitOptions);

[STL1_4state_pars,STL1_4state_likelihood] = ...
 STL1_4state_MLE.maximizeLikelihood(STL1_4state_pars,fitOptions);

% Update parameters:
for j=1:length(Model_pars)
    Model_MLE.parameters{j,2} = Model_pars(j);
end

for k=1:length(STL1_pars)
    STL1_MLE.parameters{k,2} = STL1_pars(k);
end

for l=1:length(STL1_4state_pars)
    STL1_4state_MLE.parameters{l,2} = STL1_4state_pars(l);
end

% Make plots of the parameter fits from the MLEs:
%Model_MLE.makeFitPlot
STL1_MLE.makeFitPlot
STL1_4state_MLE.makeFitPlot

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's tinker with the starting parameters of both models and try again:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of our models for refitting:
Model_MLE_refit = Model_MLE;
STL1_MLE_refit = STL1_MLE;
STL1_MLE_refit_4state = STL1_MLE_4state;

% Adjust the starting parameters:
%Model_MLE_refit.parameters = ({'kon',0.1;'koff',1;'kr',5;'gr',0.5});
STL1_MLE_refit.parameters = ({'koff',1;'kr',5;'gr',0.5; ...
                              'a0',0.1;'a1',0.01;'r1',5;'r2',1e-31});
STL1_MLE_refit_4state.parameters = ({'k12',0.5; 'k23',0.5; 'k34',0.5;...
                      'k21',0.001; 'k32',0.0001; 'k43',5e-08; ...
                      'kr1',10e-11; 'kr2',0.1; 'kr3',10; 'kr4',0.1; ...
                      'dr',0.5; 'a0',10; 'a1',5e-06; 'r1',10; 'r2',5e-16});

Model_MLE_refit.tSpan = 0:5:60;
STL1_MLE_refit.tSpan = 0:5:60;
STL1_MLE_refit_4state.tSpan = 0:5:60;

% Have SSIT approximate the steady state for initial distribution:
Model_MLE_refit.fspOptions.initApproxSS =true;
Model_MLE_refit = ...
    Model_MLE_refit.formPropensitiesGeneral('Model_MLE_refit',true);

STL1_MLE_refit.fspOptions.initApproxSS =true;
STL1_MLE_refit = ...
    STL1_MLE_refit.formPropensitiesGeneral('STL1_MLE_refit',true);

STL1_MLE_refit_4state.fspOptions.initApproxSS =true;
STL1_MLE_refit_4state = ...
STL1_MLE_refit_4state.formPropensitiesGeneral('STL1_MLE_refit_4state',true);

% Solve the FSP analysis:
[Model_MLE_refit_FSPsoln,Model_MLE_refit.fspOptions.bounds] = ...
    Model_MLE_refit.solve;   

[STL1_MLE_refit_FSPsoln,STL1_MLE_refit.fspOptions.bounds] = ...
    STL1_MLE_refit.solve;  

[STL1_MLE_refit_4state_FSPsoln,STL1_MLE_refit_4state.fspOptions.bounds] = ...
    STL1_MLE_refit_4state.solve; 

% Format new parameters:
Modelpars_refit = cell2mat(Model_MLE_refit.parameters(1:4,2));
STL1pars_refit = cell2mat(STL1_MLE_refit.parameters(1:7,2));
STL1pars_refit_4state = cell2mat(STL1_MLE_refit_4state.parameters(1:15,2));

% Maximize the likelihood:
[Modelpars_refit,Model_likelihood_refit] = ...
    Model_MLE_refit.maximizeLikelihood(Modelpars_refit, fitOptions);

[STL1pars_refit,STL1_likelihood_refit] = ...
    STL1_MLE_refit.maximizeLikelihood(STL1pars_refit, fitOptions);

[STL1pars_refit_4state,STL1_likelihood_refit_4state] = ...
STL1_MLE_refit_4state.maximizeLikelihood(STL1pars_refit_4state,fitOptions);

% Update the parameters:
for j=1:length(Modelpars_refit)
    Model_MLE_refit.parameters{j,2} = Modelpars_refit(j);
end
for k=1:length(STL1pars_refit)
    STL1_MLE_refit.parameters{k,2} = STL1pars_refit(k);
end
for l=1:length(STL1pars_refit_4state)
    STL1_MLE_refit_4state.parameters{l,2} = STL1pars_refit_4state(l);
end

% Make plots of the new model parameter fits from the MLEs:
Model_MLE_refit.makeFitPlot
STL1_MLE_refit.makeFitPlot
STL1_MLE_refit_4state.makeFitPlot
