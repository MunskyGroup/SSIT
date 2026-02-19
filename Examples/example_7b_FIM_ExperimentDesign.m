%% SSIT/Examples/example_7b_FIM_ExperimentDesign

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2.9: FIM Optimality Criteria and Experiment Design
%   * Use Fisher information results to design a next round of experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries:
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions  
% from example_4_SolveSSITModels_FSP, and loaded experimental data from 
% example_8_LoadingandFittingData_DataLoading

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_8_LoadingandFittingData_DataLoading

%% Load pre-run results:
% load('example_8_LoadingandFittingData_DataLoading.mat')

% Make a copy of our model:
STL1_4state_design = STL1_4state_data;

% Compute FIM results:
fimResults = STL1_4state_design.computeFIM([],'log',[]);

% Get the number of cells from loaded experimental data using 'nCells':
cellCounts_data = STL1_4state_design.dataSet.nCells*...
                 ones(size(STL1_4state_design.tSpan));

%% Experiment Design
% Find the FIM-based designs for a total cells 

% Compile and store propensities:
STL1_4state_design = ...
    STL1_4state_design.formPropensitiesGeneral('STL1_4state_design');

% Compute the optimal number of cells from the FIM results using different 
% design criteria:  `Trace' maximizes the trace of the FIM; 
% `DetCovariance' minimizes the expected determinant of MLE covariance; 
% `Smallest Eigenvalue' maximizes the smallest e.val of the FIM; and 
% `TR[$<i_1>,<i_2>$,...]' maximizes the determinant of the FIM for the 
% specified indices.  The latter is shown for different parameter 
% combinations, where `Tr[9:10]' are the mRNA-specific parameters `dr' and 
% `kr' (degradation and transcription, respectively).  All other parameters 
% are assumed to be known and fixed.
nCol = sum(cellCounts);
nTotal = nCol(1);
nCellsOpt_detCov = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'DetCovariance');
nCellsOpt_trace = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'Trace');
nCellsOpt_tr = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'tr[1:10]');
nCellsOpt_tr1 = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'tr[11:15]');
nCellsOpt_trR = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'tr[9:10]');

%% Make a bar chart to compare the different designs
% Find which x positions correspond to time=30 and time=60 for off-setting:
t = STL1_4state_design.tSpan;
x = 1:size(t,2);
idx = ismember(t, [30 60]);        

% Build custom x-locations for series that need separation
x_m02 = x;  x_m03(idx) = x_m02(idx) - 0.5;
x_p02 = x;  x_p02(idx) = x_p02(idx) + 0.5;

f = figure;
bar(x,      nCellsOpt_trace,  0.5); hold on
bar(x_m02,  nCellsOpt_trR,    0.5);
bar(x_p02,  nCellsOpt_tr1,    0.5);
bar(x,      nCellsOpt_detCov, 0.5);
bar(x_m02,  nCellsOpt_tr,     0.5);

set(gca,'XTick',x,'XTickLabel',t,'FontSize',16)
title('4-state STL1 (FIM Optimal Designs)','FontSize',24)
xlabel('Time (min)','FontSize',20)
ylabel('Number of cells','FontSize',20)
legend('Trace,Tr[1] Designs','Tr[9:10] Design','Tr[11:15] Design', ...
       'DetCov,\lambda Designs','Tr[1:10] Design','Location','northeast')