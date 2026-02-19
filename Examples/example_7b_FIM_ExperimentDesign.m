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

% Compile and store propensities:
STL1_4state_design = ...
    STL1_4state_design.formPropensitiesGeneral('STL1_4state_design');

%% Experiment Design
% Find the FIM-based designs for a total cells 

% Compute the optimal number of cells from the FIM results using different 
% design criteria:  `Trace' maximizes the trace of the FIM; 
% `DetCovariance' minimizes the expected determinant of MLE covariance; 
% `Smallest Eigenvalue' maximizes the smallest e.val of the FIM; and 
% `TR[$<i_1>,<i_2>$,...]' maximizes the determinant of the FIM for the 
% specified indices.  The latter is shown for different parameter 
% combinations, where `Tr[9:13]' are the mRNA-specific parameters `dr' and 
% `kr1',`kr2',`kr3', and `kr4' (degradation and transcription reactions).  
% All other parameters are assumed to be known and fixed.
nCol = sum(cellCounts);
nTotal = nCol(1);
nCellsOpt_detCov = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'DetCovariance');
nCellsOpt_trace = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'Trace');
nCellsOpt_tr = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'tr[1:8]');
nCellsOpt_trR = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'tr[9:13]');
nCellsOpt_trI = ...
    STL1_4state_design.optimizeCellCounts(fimResults,nTotal,'tr[14:18]');

%% Make a bar chart to compare the different designs
% Find which x positions correspond to time=30 and time=60 for off-setting:
t = STL1_4state_design.tSpan;
x = 1:size(t,2);   

f = figure;
bar(x,  nCellsOpt_trace,   0.4); hold on
bar(x,  nCellsOpt_detCov,  0.4);
bar(x,  nCellsOpt_trI,     0.4);
bar(x+0.2,  nCellsOpt_trR, 0.4);
bar(x-0.2,  nCellsOpt_tr,  0.4);

set(gca,'XTick',x,'XTickLabel',t,'FontSize',16)
title('4-state STL1 (FIM Optimal Designs)','FontSize',24)
xlabel('Time (min)','FontSize',20)
ylabel('Number of cells','FontSize',20)
legend('Trace Design','DetCov Design', 'Tr[14:18] Design',...
        'Tr[9:13] Design', 'Tr[1:8] Design', 'Location', 'northeast')