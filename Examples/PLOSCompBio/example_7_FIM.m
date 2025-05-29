%% example_7_FIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.3: Sensitivity analysis and Fisher Information Matrix (Part II)
% Example script to set up and solve the FSP-FIM matrix with  
% partial observations and probabilistic distortion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Load our models described in example_1_CreateSSITModels, computed FSP 
% solutions using example_4_SolveSSITModels_FSP, and computed sensitivities 
% from example_6_SensitivityAnalysis

% clear
% close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_6_SensitivityAnalysis

% View model summaries
Model_FSP.summarizeModel
STL1_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Compute the Fisher Information Matrix for the bursting gene model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the bursting gene model with solved sensitivities:
Model_FIM = Model_sens;

%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
Model_FIM = Model_sens;
Model_fimResults = Model_FIM.computeFIM(Model_sensSoln.sens); 

% Generate a count of measured cells (in place of real data):
Model_cellCounts = 10*ones(size(Model_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[Model_fimTotal,Model_mleCovEstimate,Model_fimMetrics] = ...
    Model_FIM.evaluateExperiment(Model_fimResults,Model_cellCounts)

% Plot the FIMs:
fig7 = figure(7);clf; set(fig7,'Name',...
    'Fim-Predicted Uncertainty Ellipses');
Model_FIM.plotMHResults([],Model_fimTotal,'log',[],fig7)
legend('FIM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Compute the Fisher Information Matrix for the STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the time-varying STL1 yeast model with solved 
% sensitivities:
STL1_FIM = STL1_sens;

%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
STL1_FIM = STL1_sens;
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); 

% Generate a count of measured cells (in place of real data):
STL1_cellCounts = 10*ones(size(STL1_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

% Plot the FIMs:
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