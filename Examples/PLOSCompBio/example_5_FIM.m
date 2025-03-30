%% example_5_FIMCalculation
% Example script to set up and solve the FSP-FIM matrix  
% with partial observations and probabilistic distortion.
clear
close all
addpath(genpath('../src'));

%% Preliminaries
% Load our models described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
example_2_SolveSSITModels_FSP
Model.summarizeModel
STL1Model.summarizeModel

% Make copies of Model and STL1 Model
Model_sens = Model;
STL1_sens = STL1Model;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity
Model_sens.solutionScheme = 'fspSens'; 
STL1_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem
[sensSoln,bounds] = Model_sens.solve(FSPsoln.stateSpace); 
[STL1_sensSoln,STL1_bounds] = STL1_sens.solve(STL1_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis
% Model:
fig1 = figure(1);clf; set(fig1,'Name','Marginal Sensitivity, offGene');
fig2 = figure(2);clf; set(fig2,'Name','Marginal Sensitivity, onGene');
fig3 = figure(3);clf; set(fig3,'Name','Marginal Sensitivity, mRNA');
Model_sens.makePlot(sensSoln,'marginals',[],false,...
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
fimResults = Model_FIM.computeFIM(sensSoln.sens); 

% Generate a count of measured cells (in place of real data)
cellCounts = 10*ones(size(Model_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[fimTotal,mleCovEstimate,fimMetrics] = ...
    Model_FIM.evaluateExperiment(fimResults,cellCounts)

% Plot the FIMs
fig7 = figure(7);clf; set(fig7,'Name',...
    'Fim-Predicted Uncertainty Ellipses');
Model_FIM.plotMHResults([],fimTotal,'log',[],fig7)
legend('FIM')

%% STL1 Model:
% Compute the FIM
STL1_FIM = STL1_sens;
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); 

% Generate a count of measured cells (i.e., in the case of real data,  
% the number of cells measured in the experiment)
STL1_cellCounts = 10*ones(size(STL1_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

% Plot the FIMs
fig8 = figure(5);clf; set(fig8,'Name',...
     'Fim-Predicted Uncertainty Ellipses');
STL1_FIM.plotMHResults([],STL1_fimTotal,'log',[],fig8)
legend('FIM')

% Notice if you use change the 'koff' parameter to 0.2 (where 'kon' is also 
% 0.2) and recompute the FSP solutions, then something happens to the STIL1
% Model. See the difference in fimMetrics between the basic bursting gene
% Model and the STL1 Model becomes:

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

% The determinant indicates a lack of identifiability between parameters 
% in the STL1 Model, which unlike the base model has a time-varying input 
% signal. The 'kon' and 'koff' parameters are strongly correlated and for 
% these values, would require directly observed gene data information for 
% FIM to avoid rank deficiency.
