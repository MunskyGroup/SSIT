%% example_15_PDO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Use a probability distribution operator (PDO) to handle distortion 
%     of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions  
% from example_4_SolveSSITModels_FSP, sensitivities computed in 
% example_6_SensitivityAnalysis, FIM results from example_7_FIM,  
% loaded data from example_8_LoadingandFittingData_DataLoading, and
% Metropolis-Hastings results from example_10_LoadingandFittingData_MHA
%clear
%close all
addpath(genpath('../../'));
addpath(genpath('tmpPropensityFunctions'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_6_SensitivityAnalysis
% example_7_FIM
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE
% example_10_LoadingandFittingData_MHA

% View model summary:
STL1_4state_MH.summarizeModel

% Create a copy of the STL1 model for PDO:
STL1_4state_PDO = STL1_4state_MH_FIM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Apply an affine Poisson conditional probability distribution 
%     operator (PDO) to transform parameter probabilities computed   
%     from average intensity data according to parameter probabilities 
%     computed from mRNA spot count data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figNew = figure;
plotColors = struct('scatter', [0.2, 0.6, 1], ...
           'ellipseFIM', 'r', ...
           'ellipseMH', 'g--', ...
           'marker', [0.1, 0.1, 0.1]);

sig_log10 = 2*ones(1,15);

fimTotal = ...
    STL1_4state_PDO.evaluateExperiment(fimResults,...
    STL1_4state_PDO.dataSet.nCells,diag(sig_log10.^2));

STL1_4state_PDO.plotMHResults(STL1_FIM_MHResults,[fimTotal],'log',[],...
                              figNew,plotColors)

%%

% Find and store the total number of cells in your data set (already
% computed by SSIT when data was loaded in example_8:
nTotal = sum(STL1_4state_PDO.dataSet.nCells);

% Compute the optimal number of cells from the FIM results computed in 
% example_7_FIM using the min. inv determinant <x^{-1}> 
% (all other parameters are known and fixed)
nCellsOpt = ...
    STL1_4state_PDO.optimizeCellCounts(fimResults,...
                                       nTotal,'tr[1:15]');

nCellsOptAvail = min(nCellsOpt,STL1_4state_PDO.dataSet.nCells)
fimOpt = STL1_4state_PDO.evaluateExperiment(fimResults,nCellsOpt,...
                                            diag(sig_log10.^2));
fimOptAvail = STL1_4state_PDO.evaluateExperiment(fimResults,...
                                        nCellsOptAvail,diag(sig_log10.^2));
figOpt = figure;
STL1_4state_PDO.plotMHResults(STL1_FIM_MHResults,[fimOpt,fimTotal],...
                              'log',[],figOpt,plotColors);
figOptAvail = figure;
STL1_4state_PDO.plotMHResults(STL1_FIM_MHResults,[fimOptAvail,fimTotal],...
                              'log',[],figOptAvail,plotColors);

f = figure;
set(f,'Position',[616   748   412   170])
bar([1:16],STL1_4state_PDO.dataSet.nCells,0.45)
hold on
bar([1:16]+0.5,nCellsOpt,0.45)
set(gca,'xtick',[1:16]+0.25,'xticklabel',STL1_4state_PDO.dataSet.times,...
                                         'fontsize',16,'ylim',[0,7000])
legend('Intuitive Design','Optimal Design')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PDO Calculations
%% Ex(1): Calibrate PDO from nuclear mRNA count data
% Calibrate the PDO from empirical data. Here, the number of spots has
% been measured using different assays in data columns 'nTotal' for the
% 'true' data set and in the columns 'nSpots0' for a different label or
% 'intens1' for the integrated intensity.  We calibrate two different 
% PDOs for this case. In both cases, we assume an 'AffinePoiss' PDO  
% where the obervation probability is a Poisson distribution where the   
% mean value is affine linearly related to the true value: 
% P(y|x) = Poiss(a0 + a1*x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STL1_4state_PDO = ...
  STL1_4state_PDO.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
  {'mRNA'},{'RNA_STL1_total_TS3Full'},{'RNA_STL1_nuc_TS3Full'},...
   'AffinePoiss',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Calibrate PDO from average intensity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STL1_4state_PDO_intens = STL1_4state_PDO;
STL1_4state_PDO_intens = ...
STL1_4state_PDO_intens.calibratePDO(...
        'data/filtered_data_2M_NaCl_Step.csv',...
        {'mRNA'},{'RNA_STL1_total_TS3Full'},{'STL1_avg_int_TS3Full'},...
        'AffinePoiss',true,[1,4000,80]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM + PDO analyses
%   * Analyze FIM with PDO for nuclear mRNA count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fimsPDO = STL1_4state_PDO.computeFIM([],'log');
fimPDO = STL1_4state_PDO.evaluateExperiment(fimsPDO,nCellsOpt,...
                                            diag(sig_log10.^2));

nCellsOptPDO = STL1_4state_PDO.optimizeCellCounts(fimsPDO,nTotal,...
                                                  'tr[1:15]');

figPDO = figure;
STL1_4state_PDO.plotMHResults(STL1_FIM_MHResults,...
                             [fimPDO,fimTotal,fimOpt],'log',[],figPDO);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM + PDO analyses
%   * Analyze FIM with PDO for average intensity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fimsPDOintens = STL1_4state_PDO_intens.computeFIM([],'log');
fimPDOintens = STL1_4state_PDO_intens.evaluateExperiment(fimsPDOintens,...
                                            nCellsOpt,diag(sig_log10.^2));

nCellsOptPDOintens = STL1_4state_PDO_intens.optimizeCellCounts(...
                            fimsPDOintens,nTotal,'tr[1:15]');

figintens = figure;
STL1_4state_PDO_intens.plotMHResults(STL1_FIM_MHResults,...
                        [fimPDOintens,fimTotal,fimOpt],'log',[],figintens);

%% Save PDO models + results:
saveNames = unique({'STL1_4state_PDO'
    'fimTotal'
    'nTotal'
    'nCellsOpt'
    'nCellsOptAvail'
    'fimOpt'
    'fimOptAvail'
    'fimsPDO'
    'fimPDO'
    'nCellsOptPDO'
    'STL1_4state_PDO_intens'
    'fimsPDOintens'
    'fimPDOintens'
    'nCellsOptPDOintens'
    });
    
save('example_15_PDO',saveNames{:})
