%% example_6_LoadingandFittingData_MLE
% Example script to demonstrate how to load and fit
% experimental data using maximum likelihood estimates (MLEs)
addpath('../../src');

%% Preliminaries
% Load our STL1 Model described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
%% Comment out the following 3 lines if example_2 has already been run:
clear
close all
example_2_SolveSSITModels_FSP

% View model summaries
Model_FSP.summarizeModel
STL1_FSP.summarizeModel

% Make new copies of our models
ModelReal = Model_FSP;
STL1Real = STL1_FSP;

% Set fitOptions, with the maximum allowable number of iterations to fit
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them)
% and convert from cell to double
Modelpars = cell2mat(Model.parameters(1:4,2));
STL1pars = cell2mat(STL1.parameters(1:7,2));

%% Load experimental data, matching the species name of the model ('mRNA') 
% to the appropriate column in the data file ('rna')
ModelReal = ModelReal.loadData('data/STL1.csv',{'mRNA','rna'});
STL1Real = STL1Real.loadData('data/STL1.csv',{'mRNA','rna'});

% These plots are unnecessary, as the model parameters have not been fit to
% the data yet.  However, it illustrates the improvement due to MLEs later:
%ModelReal.makeFitPlot
%STL1Real.makeFitPlot

%% Compute the MLE
% ModelReal:
[Modelpars,Model_likelihood] = ModelReal.maximizeLikelihood(Modelpars,...
                                                            fitOptions);
% STL1Real:
[STL1pars,STL1_likelihood] = ...
 STL1Real.maximizeLikelihood(STL1pars,fitOptions);

% Update ModelReal parameters:
for i=1:length(Modelpars)
    ModelReal.parameters{i,2} = Modelpars(i);
end
% Update STL1Real parameters
for j=1:length(STL1pars)
    STL1Real.parameters{j,2} = STL1pars(j);
end

% Make plots of the model parameter fits from the MLEs
ModelReal.makeFitPlot
STL1Real.makeFitPlot

%% Let's tinker with the starting parameters of the STL1 Model 
%% and try again:
STL1Real_refit = STL1Real;

STL1Real_refit.parameters = ({'koff',30;'kr',100;'gr',0.005; ...
                         'a0',0.01;'a1',1;'r1',0.4;'r2',0.1});

STL1Real_refit.tSpan = 0:5:60;

% Have SSIT approximate the steady state for initial distribution
STL1Real_refit.fspOptions.initApproxSS =true;
STL1Real_refit = ...
    STL1Real_refit.formPropensitiesGeneral('STL1Real_refit',true);

% Solve the FSP analysis
[STL1Real_refit_FSPsoln,STL1Real_refit.fspOptions.bounds] = ...
    STL1Real_refit.solve;  

% Format new parameters
STL1pars_refit = cell2mat(STL1Real_refit.parameters(1:7,2));

% Maximize the likelihood
[STL1pars_refit,STL1_likelihood1] = ...
    STL1Real_refit.maximizeLikelihood(STL1pars_refit, fitOptions);

% Update STL1Real parameters
for j=1:length(STL1pars_refit)
    STL1Real_refit.parameters{j,2} = STL1pars_refit(j);
end

% Make plots of the new model parameter fits from the MLEs
STL1Real_refit.makeFitPlot