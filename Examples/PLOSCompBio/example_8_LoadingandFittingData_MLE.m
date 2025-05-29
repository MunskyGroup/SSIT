%% example_8_LoadingandFittingData_MLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Loading and fitting time-varying STL1 yeast data 
%%     * Data loading and handling 
%%     * Maximize the likelihood L(D|theta) and use the maximum 
%%       likelihood estimate (MLE) to fit the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels and  
% computed FSP solutions from example_4_SolveSSITModels_FSP
%clear
%close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP

% View model summary:
STL1_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load time-varying STL1 yeast data and filter by experimental 
%% condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a new copy of our model:
STL1_data = STL1_FSP;

%% Load the experimental data, matching the species name of the model  
%% ('mRNA') to the appropriate column in the data file ('rna') and filter 
%% it so that only data from experiment 1, replica 1 is loaded:

% STL1_data = STL1_data.loadData('data/Hog1_exp_rep_CY5_combined.csv',...
%                               {'mRNA','rna'},{'experiment',1;'replica',1});
STL1_data = STL1_data.loadData('data/STL1.csv',...
                               {'mRNA','rna'});

% This plot is unnecessary, as the model parameters have not been fit to
% the data yet.  However, it illustrates the improvement to come later:
STL1_data.makeFitPlot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit experimental data using maximum likelihood estimates (MLEs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of STL1_data for MLE:
STL1_MLE = STL1_data;

% Set fitOptions, with the maximum allowable number of iterations to fit
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
STL1pars = cell2mat(STL1_MLE.parameters(1:7,2));

%% Compute the MLE
[STL1pars,STL1_likelihood] = ...
 STL1_MLE.maximizeLikelihood(STL1pars,fitOptions);

% Update STL1 parameters
for j=1:length(STL1pars)
    STL1_MLE.parameters{j,2} = STL1pars(j);
end

% Make plots of the model parameter fits from the MLEs
STL1_MLE.makeFitPlot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's tinker with the starting parameters of the STL1 Model and try 
%% again:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of STL1_MLE for refitting:
STL1_MLE_refit = STL1_MLE;

STL1_MLE_refit.parameters = ({'koff',30;'kr',100;'gr',0.005; ...
                              'a0',0.01;'a1',1;'r1',0.4;'r2',0.1});

STL1_MLE_refit.tSpan = 0:5:60;

% Have SSIT approximate the steady state for initial distribution
STL1_MLE_refit.fspOptions.initApproxSS =true;
STL1_MLE_refit = ...
    STL1_MLE_refit.formPropensitiesGeneral('STL1_MLE_refit',true);

% Solve the FSP analysis
[STL1_MLE_refit_FSPsoln,STL1_MLE_refit.fspOptions.bounds] = ...
    STL1_MLE_refit.solve;  

% Format new parameters
STL1pars_refit = cell2mat(STL1_MLE_refit.parameters(1:7,2));

% Maximize the likelihood
[STL1pars_refit,STL1_likelihood_refit] = ...
    STL1_MLE_refit.maximizeLikelihood(STL1pars_refit, fitOptions);

% Update STL1Real parameters
for j=1:length(STL1pars_refit)
    STL1_MLE_refit.parameters{j,2} = STL1pars_refit(j);
end

% Make plots of the new model parameter fits from the MLEs
STL1_MLE_refit.makeFitPlot