%% SSIT/Examples/example_7_FIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Fisher Information Matrix (FIM)
%   * Set up and solve the FSP-FIM matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels and computed FSP solutions 
% from example_4_SolveSSITModels_FSP

% clear
% close all

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP

%% Load pre-computed FSP solutions:
load('example_4_SolveSSITModels_FSP.mat')

% View model summaries:
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Compute the Fisher Information Matrix for the bursting gene model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
Model_fimResults = Model.computeFIM(scale='log'); 

% Generate a count of measured cells (in place of real data):
Model_cellCounts = 100*ones(size(Model.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[Model_fimTotal,Model_mleCovEstimate,Model_fimMetrics] = ...
    Model.evaluateExperiment(Model_fimResults,Model_cellCounts)

theta0 = [Model.parameters{:,2}];

Model.plotFIMResults(Model_fimTotal, 'log', Model.parameters,...
                         theta0, PlotEllipses=true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Compute the Fisher Information Matrix for the STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute FIMs using FSP sensitivity results
% Compute the FIM:
STL1_fimResults = STL1.computeFIM(scale='log'); 

% Generate a count of measured cells (in place of real data):
STL1_cellCounts = 100*ones(size(STL1.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

STL1_theta0 = [STL1.parameters{:,2}];

% Plot the FIMs:
STL1.plotFIMResults(STL1_fimTotal, 'log', STL1.parameters,...
                        STL1_theta0, PlotEllipses=true,...
                        EllipsePairs=[1 3; 1 4; 3 4; 7 8; 3 6; 2 7]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Compute the FIM for the 4-state STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define indices of free parameters for FIM sub matrix. (Hog1 input signal
%% params are experimentally known (thus fixed) and all others are free:
freePars = 1:13;

% Specify time points and numbers of cells for experiments
STL1_4state.tSpan = [0,1,2,4,6,8,10:5:55];
cellCounts = 1000*ones(1,16);

%% Compute FIMs using FSP sensitivity results
% Compute the full FIM:
fims_full = STL1_4state.computeFIM(scale='log',observed='mRNA'); % full FIM

% Compute the FIM sub matrix for free parameters:
fims_free = STL1_4state.computeFIM(scale='log',observed='mRNA',...
    freePars=freePars); % FIM sub matrix

% - Or, get the number of cells using 'nCells':
% STL1_4state_cellCounts = ...
% STL1_4state_data.dataSet.nCells*ones(size(STL1_4state.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set):
[fimTotal_full,mleCovEstimate_full,fimMetrics_full] = ...
    STL1_4state.evaluateExperiment(fims_full,cellCounts)

[fimTotal_free,mleCovEstimate_free,fimMetrics_free] = ...
    STL1_4state.evaluateExperiment(fims_free,cellCounts)

%% Plot the FIMs (full):
f1 = figure(11);
f2 = figure(12);
STL1_4state.plotFIMResults(fimTotal_full,'log',STL1_4state.parameters,...
    PlotEllipses=true, EllipseFigure=f1,...
    EllipsePairs=[1 6; 2 3; 4 5; 6 13], FigureHandle=f2,...
    Colors=struct('EllipseColors',[0.2 0.6 0.9],...
    'CenterSquare',[0.96,0.47,0.16]));

% Plot the FIMs (free):
f3 = figure(13);
STL1_4state.plotFIMResults(fimTotal_free, 'log',...
    STL1_4state.parameters(1:13),PlotEllipses=true,EllipseFigure=f1,...
    EllipsePairs=[1 6; 2 3; 4 5; 6 13],FigureHandle=f3,...
    Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]));

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
%          det: 4.8027e+24
%        trace: 1.6623e+07
%    minEigVal: 1.6309e+05

% STL1_fimMetrics = 
%          det: 9.9840e+00
%        trace: 3.5112e+04
%    minEigVal: 2.5775e-11

%% Full FIM
% STL1_4state_fimMetrics = 
%          det: 1.7388e+09
%        trace: 2.4234e+05
%    minEigVal: -3.4038e-14

%% FIM sub matrix (free parameters):
% STL1_4state_fimMetrics_free =  
%           det: 2.7363e+29
%         trace: 7.2778e+04
%     minEigVal: 2.9656e-01


%% Save models & FIM results
saveNames = unique({'Model'
    'Model_fimResults'
    'Model_cellCounts'
    'Model_fimTotal'
    'Model_mleCovEstimate'
    'Model_fimMetrics'
    'STL1'
    'STL1_fimResults'
    'STL1_cellCounts'
    'STL1_fimTotal'
    'STL1_mleCovEstimate'
    'STL1_fimMetrics'
    'STL1_4state'
    'STL1_4state_fimResults_full'
    'STL1_4state_fimResults_free'
    'cellCounts'
    'STL1_4state_fimTotal_full'
    'STL1_4state_fimTotal_free'
    'STL1_4state_mleCovEstimate_full'
    'STL1_4state_mleCovEstimate_free'
    'STL1_4state_fimMetrics_full'
    'STL1_4state_fimMetrics_free'
    });
    
save('example_7_FIM',saveNames{:})
