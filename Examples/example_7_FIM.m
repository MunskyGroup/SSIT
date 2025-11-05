%% example_7_FIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.3: Sensitivity analysis and Fisher Information Matrix 
%% (Part II)
%   * Set up and solve the FSP-FIM matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels, computed FSP solutions 
% from example_4_SolveSSITModels_FSP, and computed sensitivities from 
% example_6_SensitivityAnalysis

% clear
% close all
addpath(genpath('../'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_6_SensitivityAnalysis

%% Load pre-computed sensitivities:
load('example_6_SensitivityAnalysis.mat')

% View model summaries:
Model_sens.summarizeModel
STL1_sens.summarizeModel
STL1_4state_sens.summarizeModel

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
fig1 = figure(12);clf; set(fig1,'Name',...
    'Fim-Predicted Uncertainty Ellipses');
Model_FIM.plotMHResults([],Model_fimTotal,'log',[],fig1)
legend('FIM')

Model_FIM.plotFIMResults(Model_fimTotal,Model_FIM.parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Compute the Fisher Information Matrix for the STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the time-varying STL1 yeast model with solved 
% sensitivities:
STL1_FIM = STL1_sens;

%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); 

% Generate a count of measured cells (in place of real data):
STL1_cellCounts = 10*ones(size(STL1_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

% Plot the FIMs:
fig2 = figure(13);clf; set(fig2,'Name',...
     'Fim-Predicted Uncertainty Ellipses');
STL1_FIM.plotMHResults([],STL1_fimTotal,'log',[],fig2)
legend('FIM')

STL1_FIM.plotFIMResults(STL1_fimTotal,STL1_FIM.parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Compute the FIM for the 4-state STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the 4-state time-varying STL1 yeast model with solved 
% sensitivities:
STL1_4state_FIM = STL1_4state_sens;

%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
STL1_4state_fimResults = ...
    STL1_4state_FIM.computeFIM(); 

% Generate a count of measured cells - or get the number of cells using 
% 'nCells' - e.g., "STL1_4state_data.dataSet.nCells", which becomes:
% STL1_4state_cellCounts = ...
% STL1_4state_data.dataSet.nCells*ones(size(STL1_4state_FIM.tSpan));
STL1_4state_cellCounts = 10*ones(size(STL1_4state_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[STL1_4state_fimTotal,STL1_4state_mleCovEstimate,...
    STL1_4state_fimMetrics] = ...
    STL1_4state_FIM.evaluateExperiment(STL1_4state_fimResults,...
                                       STL1_4state_cellCounts)

% Plot the FIMs:
fig3 = figure(14);clf; set(fig3,'Name',...
     'Fim-Predicted Uncertainty Ellipses');
STL1_4state_FIM.plotMHResults([],STL1_4state_fimTotal,'log',[],fig3)
legend('FIM')

STL1_4state_FIM.plotFIMResults(STL1_4state_fimTotal,...
    STL1_4state_FIM.parameters)

%%
% Note:  If detI(θ)=0, then at least one eigenvalue 𝜆𝑘=0. That means the 
% FIM is rank-deficient, so there is at least one non-trivial linear 
% combination of parameters whose variance (via the Cramér–Rao bound) is 
% infinite. That direction is locally non-identifiable at 𝜃; the likelihood 
% is flat in that direction.  If detI(θ) is nonzero but extremely small, 
% that usually means one or more eigenvalues are tiny (but not exactly 
% zero). Then there is practical non-identifiability or very strong 
% parameter correlation: those directions in parameter space are only very 
% weakly constrained by your experiment.

% Model_fimMetrics = 
% 
%   struct with fields:
% 
%           det: 3.978068235240957e+08
%         trace: 2.291624663469069e+04
%     minEigVal: 1.859230327938760e-01

% STL1_fimMetrics = 
% 
%   struct with fields:
% 
%           det: 3.956395658594882e-14
%         trace: 3.491622747290549e+04
%     minEigVal: -7.945218142991907e-14

% STL1_4state_fimMetrics = 
% 
%   struct with fields:
% 
%           det: 1.585910799906740e+24
%         trace: 1.276663861682894e+11
%     minEigVal: 2.408883953969782e-11


%% Save models & FIM results
saveNames = unique({'Model_FIM'
    'Model_fimResults'
    'Model_cellCounts'
    'Model_fimTotal'
    'Model_mleCovEstimate'
    'Model_fimMetrics'
    'STL1_FIM'
    'STL1_fimResults'
    'STL1_cellCounts'
    'STL1_fimTotal'
    'STL1_mleCovEstimate'
    'STL1_fimMetrics'
    'STL1_4state_FIM'
    'STL1_4state_fimResults'
    'STL1_4state_cellCounts'
    'STL1_4state_fimTotal'
    'STL1_4state_mleCovEstimate'
    'STL1_4state_fimMetrics'
    });
    
save('example_7_FIM',saveNames{:})