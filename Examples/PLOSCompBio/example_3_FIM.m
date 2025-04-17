%% example_5_FIMCalculation
% Example script to set up and solve the FSP-FIM matrix  
% with partial observations and probabilistic distortion.
addpath(genpath('../../src'));

%% Preliminaries
% Load our models described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
%% Comment out the following 3 lines if example_4_FIM has already been run:
clear
close all
example_3_LoadingandFittingData_MLE

% View model summaries
ModelReal.summarizeModel
STL1Real_refit.summarizeModel

Model_sens = ModelReal;
STL1_sens = STL1Real_refit;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity
Model_sens.solutionScheme = 'fspSens'; 
STL1_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem
[Model_sensSoln,Model_bounds] = Model_sens.solve(Model_FSPsoln.stateSpace); 
[STL1_sensSoln,STL1_bounds] = STL1_sens.solve(STL1_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis
% Model:
fig1 = figure(1);clf; set(fig1,'Name','Marginal Sensitivity, offGene');
fig2 = figure(2);clf; set(fig2,'Name','Marginal Sensitivity, onGene');
fig3 = figure(3);clf; set(fig3,'Name','Marginal Sensitivity, mRNA');
Model_sens.makePlot(Model_sensSoln,'marginals',[],false,...
                    [fig1,fig2,fig3],{'b','linewidth',2})
% STL1 Model:
fig4 = figure(4);clf; set(fig4,'Name','Marginal Sensitivity, offGene');
fig5 = figure(5);clf; set(fig5,'Name','Marginal Sensitivity, onGene');
fig6 = figure(6);clf; set(fig6,'Name','Marginal Sensitivity, mRNA');
STL1_sens.makePlot(STL1_sensSoln,'marginals',[],false,...
                   [fig4,fig5,fig6],{'b','linewidth',2})

%% Compute FIMs using FSP sensitivity results
%% Model:
% Compute the FIM
Model_FIM = Model_sens;
Model_fimResults = Model_FIM.computeFIM(Model_sensSoln.sens); 

% Generate a count of measured cells (in place of real data)
Model_cellCounts = 10*ones(size(Model_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[Model_fimTotal,Model_mleCovEstimate,Model_fimMetrics] = ...
    Model_FIM.evaluateExperiment(Model_fimResults,Model_cellCounts)

% Plot the FIMs
fig7 = figure(7);clf; set(fig7,'Name',...
    'Fim-Predicted Uncertainty Ellipses');
Model_FIM.plotMHResults([],Model_fimTotal,'log',[],fig7)
legend('FIM')

%% STL1 Model:
% Compute the FIM
STL1_FIM = STL1_sens;
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); 

% The number of cells measured in the STL1 experiment
STL1_cellCounts = STL1_FIM.dataSet.nCells;

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

% Plot the FIMs
fig8 = figure(5);clf; set(fig8,'Name',...
     'Fim-Predicted Uncertainty Ellipses');
STL1_FIM.plotMHResults([],STL1_fimTotal,'log',[],fig8)
legend('FIM')

%%
% Note:  Under certain circumstances, the FIM calculation will fail.  For 
% example, without the removel of the 'kon' parameter from STL1 (which has 
% no effect on STL1 unlike the basic bursting gene Model), Model_fimMetrics 
% (from the basic bursting gene Model) and STL1_fimMetrics (from the 
% time-varying input model) become:

% fimMetrics = 
% 
%   struct with fields:
% 
%           det: 1.050283981892288e+12
%         trace: 1.875101485843730e+04
%     minEigVal: 1.010862869918453e+02

% STL1_fimMetrics = 
% 
%   struct with fields:
% 
%           det: 0
%         trace: 5.569716495312674e+04
%     minEigVal: 0

% The determinant for STL1_fimMetrics indicates a lack of identifiability 
% between parameters in the STL1 Model, which unlike the base model has a 
% time-varying input signal. There is not information for the STL1 Model 
% concerning 'kon', leading to rank deficiency for the FIM.