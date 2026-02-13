%% SSIT/Examples/example_7_FIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.3: Sensitivity analysis and Fisher Information Matrix 
%% (Part II)
%   * Set up and solve the FSP-FIM matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels and computed FSP solutions 
% from example_4_SolveSSITModels_FSP

% clear
% close all
addpath(genpath('../'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP

%% Load pre-computed sensitivities:
% load('example_4_SolveSSITModels_FSP.mat')

% View model summaries:
Model_FSP.summarizeModel
STL1_FSP.summarizeModel
STL1_4state_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Compute the Fisher Information Matrix for the bursting gene model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the bursting gene model with solved sensitivities:
Model_FIM = Model_FSP;

%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
Model_FIM = Model_FSP;
Model_fimResults = Model_FIM.computeFIM([],'log',[]); 

% Generate a count of measured cells (in place of real data):
Model_cellCounts = 100*ones(size(Model_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[Model_fimTotal,Model_mleCovEstimate,Model_fimMetrics] = ...
    Model_FIM.evaluateExperiment(Model_fimResults,Model_cellCounts)

theta0 = [Model_FIM.parameters{:,2}];

Model_FIM.plotFIMResults(Model_fimTotal, 'log', Model_FIM.parameters,...
                         theta0, PlotEllipses=true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Compute the Fisher Information Matrix for the STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the time-varying STL1 yeast model with solved 
% sensitivities:
STL1_FIM = STL1_FSP;

%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
STL1_fimResults = STL1_FIM.computeFIM([],'log',[]); 

% Generate a count of measured cells (in place of real data):
STL1_cellCounts = 100*ones(size(STL1_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

STL1_theta0 = [STL1_FIM.parameters{:,2}];

% Plot the FIMs:
STL1_FIM.plotFIMResults(STL1_fimTotal, 'log', STL1_FIM.parameters,...
                        STL1_theta0, PlotEllipses=true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Compute the FIM for the 4-state STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the 4-state time-varying STL1 yeast model with solved 
% sensitivities:
STL1_4state_FIM = STL1_4state_FSP;

% Define indices of free parameters for FIM sub matrix. Here, Hog1 input 
% parameters are experimentally known (thus fixed) and all others are free:
freePars = 1:13;

%% Compute FIMs using FSP sensitivity results
% Compute the full FIM:
STL1_4state_fimResults_full = STL1_4state_FIM.computeFIM([],'log',[]); 

% Compute the FIM sub matrix for free parameters:
STL1_4state_fimResults_free = ...
    STL1_4state_FIM.computeFIM([],'log',[],freePars);

% Generate a count of measured cells:
STL1_4state_cellCounts = 1000*ones(size(STL1_4state_FIM.tSpan));

% - Or, get the number of cells using 'nCells':
% STL1_4state_cellCounts = ...
% STL1_4state_data.dataSet.nCells*ones(size(STL1_4state_FIM.tSpan));


% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[STL1_4state_fimTotal_full,STL1_4state_mleCovEstimate_full,...
    STL1_4state_fimMetrics_full] = ...
    STL1_4state_FIM.evaluateExperiment(STL1_4state_fimResults_full,...
                                       STL1_4state_cellCounts)

[STL1_4state_fimTotal_free,STL1_4state_mleCovEstimate_free,...
    STL1_4state_fimMetrics_free] = ...
    STL1_4state_FIM.evaluateExperiment(STL1_4state_fimResults_free,...
                                       STL1_4state_cellCounts)

% Plot the FIMs (6 parameter combinations):
STL1_4state_FIM.plotFIMResults(STL1_4state_fimTotal_full, 'log',...
    STL1_4state_FIM.parameters, PlotEllipses=true,...
    EllipsePairs=[1 6; 9 17; 2 3; 9 10; 8 13; 8 9]);

% Plot the FIMs (9 parameter combinations):
STL1_4state_FIM.plotFIMResults(STL1_4state_fimTotal_free, 'log',...
    STL1_4state_FIM.parameters(1:13));

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
%          det: 3.9036e+21
%        trace: 2.2181e+07
%    minEigVal: 2.6411e+03

% STL1_fimMetrics = 
%          det: 1.5059e+00
%        trace: 4.5588e+05
%    minEigVal: 6.7265e-12

%% Full FIM
% STL1_4state_fimMetrics = 
%          det: 3.1402e+00
%        trace: 9.2817e+08
%    minEigVal: 8.4878e-13

%% FIM sub matrix (free parameters):
% STL1_4state_fimMetrics_free =  
%           det: 2.7889e+44
%         trace: 6.9367e+05
%     minEigVal: 1.1381


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
    'STL1_4state_FIM'
    'STL1_4state_fimResults_full'
    'STL1_4state_fimResults_free'
    'STL1_4state_cellCounts'
    'STL1_4state_fimTotal_full'
    'STL1_4state_fimTotal_free'
    'STL1_4state_mleCovEstimate_full'
    'STL1_4state_mleCovEstimate_free'
    'STL1_4state_fimMetrics_full'
    'STL1_4state_fimMetrics_free'
    });
    
save('example_7_FIM',saveNames{:})
