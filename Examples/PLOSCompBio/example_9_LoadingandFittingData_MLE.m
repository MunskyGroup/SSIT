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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit experimental data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make new copies of our model, a simple bursting gene model and a more 
% complex bursting gene model that has a time-varying input signal that 
% turns on the STL1 gene:
Model_MLE = Model_data;
STL1_MLE = STL1_data;
% Let's see which model better fits our data...

% Set fitOptions, with the maximum allowable number of iterations to fit:
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
Modelpars = cell2mat(Model_MLE.parameters(1:4,2));
STL1pars = cell2mat(STL1_MLE.parameters(1:7,2));

%% Compute the MLEs:
[Modelpars,Model_likelihood] = ...
 Model_MLE.maximizeLikelihood(Modelpars,fitOptions);

[STL1pars,STL1_likelihood] = ...
 STL1_MLE.maximizeLikelihood(STL1pars,fitOptions);

% Update parameters:
for j=1:length(Modelpars)
    Model_MLE.parameters{j,2} = Modelpars(j);
end
for k=1:length(STL1pars)
    STL1_MLE.parameters{k,2} = STL1pars(k);
end

% Make plots of the parameter fits from the MLEs:
Model_MLE.makeFitPlot
STL1_MLE.makeFitPlot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's tinker with the starting parameters of both models and try again:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of our models for refitting:
Model_MLE_refit = Model_MLE;
STL1_MLE_refit = STL1_MLE;

% Adjust the starting parameters:
Model_MLE_refit.parameters = ({'kon',0.1;'koff',30;'kr',100;'gr',0.005});
STL1_MLE_refit.parameters = ({'koff',30;'kr',100;'gr',0.005; ...
                              'a0',0.01;'a1',1;'r1',0.4;'r2',0.1});

Model_MLE_refit.tSpan = 0:5:60;
STL1_MLE_refit.tSpan = 0:5:60;

% Have SSIT approximate the steady state for initial distribution:
Model_MLE_refit.fspOptions.initApproxSS =true;
Model_MLE_refit = ...
    Model_MLE_refit.formPropensitiesGeneral('Model_MLE_refit',true);

STL1_MLE_refit.fspOptions.initApproxSS =true;
STL1_MLE_refit = ...
    STL1_MLE_refit.formPropensitiesGeneral('STL1_MLE_refit',true);

% Solve the FSP analysis:
[Model_MLE_refit_FSPsoln,Model_MLE_refit.fspOptions.bounds] = ...
    Model_MLE_refit.solve;  

[STL1_MLE_refit_FSPsoln,STL1_MLE_refit.fspOptions.bounds] = ...
    STL1_MLE_refit.solve;  

% Format new parameters:
Modelpars_refit = cell2mat(Model_MLE_refit.parameters(1:4,2));
STL1pars_refit = cell2mat(STL1_MLE_refit.parameters(1:7,2));

% Maximize the likelihood:
[Modelpars_refit,Model_likelihood_refit] = ...
    Model_MLE_refit.maximizeLikelihood(Modelpars_refit, fitOptions);

[STL1pars_refit,STL1_likelihood_refit] = ...
    STL1_MLE_refit.maximizeLikelihood(STL1pars_refit, fitOptions);

% Update the parameters:
for j=1:length(Modelpars_refit)
    Model_MLE_refit.parameters{j,2} = Modelpars_refit(j);
end
for k=1:length(STL1pars_refit)
    STL1_MLE_refit.parameters{k,2} = STL1pars_refit(k);
end

% Make plots of the new model parameter fits from the MLEs:
Model_MLE_refit.makeFitPlot
STL1_MLE_refit.makeFitPlot