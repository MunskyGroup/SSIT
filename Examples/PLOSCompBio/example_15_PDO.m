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

example_1_CreateSSITModels  
example_4_SolveSSITModels_FSP
example_6_SensitivityAnalysis
example_7_FIM
example_8_LoadingandFittingData_DataLoading
example_9_LoadingandFittingData_MLE
example_10_LoadingandFittingData_MHA

% View model summary:
STL1_MH.summarizeModel

% Create a copy of the STL1 model for PDO:
STL1_PDO = STL1_MH;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Apply an affine Poisson conditional probability distribution 
%     operator (PDO) to transform parameter probabilities computed   
%     from average intensity data according to parameter probabilities 
%     computed from mRNA spot count data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find and store the total number of cells in your data set (already
% computed by SSIT when data was loaded in example_8:
nTotal = sum(STL1_PDO.dataSet.nCells);

% Compute the optimal number of cells from the FIM results computed in 
% example_7_FIM using the min. inv determinant <x^{-1}> 
% (all other parameters are known and fixed)
nCellsOpt = STL1_PDO.optimizeCellCounts(STL1_fimResults,...
                                        nTotal,'tr[1:7]');

sig_log10 = 2*ones(1,7); 

nCellsOptAvail = min(nCellsOpt,STL1_PDO.dataSet.nCells)
fimOpt = STL1_PDO.evaluateExperiment(STL1_fimResults,...
                                        nCellsOpt,diag(sig_log10.^2));
fimOptAvail = STL1_PDO.evaluateExperiment(STL1_fimResults,...
                                      nCellsOptAvail,diag(sig_log10.^2));
figOpt = figure;
STL1_PDO.plotMHResults(STL1_MHResults,[fimOpt,fimTotal],'log',[],figOpt);
figOptAvail = figure;
STL1_PDO.plotMHResults(STL1_MHResults,[fimOptAvail,fimTotal],'log',...
                                                         [],figOptAvail);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[1,0,1];
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];
        CH(2).LineWidth = 3;
        CH(3).Color=[0,1,1];
        CH(3).LineWidth = 3;
    end
end

f = figure;
set(f,'Position',[616   748   412   170])
bar([1:16],STL1_PDO.dataSet.nCells,0.45)
hold on
bar([1:200]+0.5,nCellsOpt,0.45)
set(gca,'xtick',[1:16]+0.25,'xticklabel',STL1_PDO.dataSet.times,...
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
STL1_PDO = ...
  STL1_PDO.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
  {'mRNA'},{'RNA_STL1_total_TS3Full'},{'RNA_STL1_nuc_TS3Full'},...
   'AffinePoiss',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Calibrate PDO from average intensity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STL1_PDO_intens = STL1_PDO;
STL1_PDO_intens = ...
STL1_PDO_intens.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
    {'mRNA'},{'RNA_STL1_total_TS3Full'},{'STL1_avg_int_TS3Full'},...
     'AffinePoiss',true,[1,4000,80]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM + PDO analyses
%   * Analyze FIM with PDO for nuclear mRNA count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fimsPDO = STL1_PDO.computeFIM([],'log');
fimPDO = STL1_PDO.evaluateExperiment(fimsPDO,...
                                     nCellsOpt,diag(sig_log10.^2));

nCellsOptPDO = STL1_PDO.optimizeCellCounts(fimsPDO,nTotal,'tr[1:7]');


figPDO = figure;
STL1_PDO.plotMHResults(STL1_MHResults,...
    [fimPDO,fimTotal,fimOpt],'log',[],figPDO);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[0,0,0];   % MH - black
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];   % MLE - black
        CH(2).LineWidth = 3;
        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
        CH(3).LineWidth = 3;
        CH(4).Color=[0,1,1];   % fimTotal - blue 
        CH(4).LineWidth = 3;
        CH(5).Color=[1,0,1];   % fimOpt - magenta
        CH(5).LineWidth = 3;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM + PDO analyses
%   * Analyze FIM with PDO for average intensity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fimsPDOintens = STL1_PDO_intens.computeFIM([],'log');
fimPDOintens = STL1_PDO_intens.evaluateExperiment(fimsPDOintens,...
                                          nCellsOpt,diag(sig_log10.^2));

nCellsOptPDO = STL1_PDO_intens.optimizeCellCounts(fimsPDOintens,...
                                                  nTotal,'tr[1:7]');


figintens = figure;
STL1_PDO_intens.plotMHResults(STL1_MHResults,...
                     [fimPDOintens,fimTotal,fimOpt],'log',[],figintens);
for i = 1:3
    for j = i:3
        subplot(3,3,(i-1)*3+j)
        CH = get(gca,'Children');
        CH(1).Color=[0,0,0];   % MH - black
        CH(1).LineWidth = 3;
        CH(2).Color=[0,0,0];   % MLE - black
        CH(2).LineWidth = 3;
        CH(3).Color=[0,0,1];   % fimPDOSpots - cyan
        CH(3).LineWidth = 3;
        CH(4).Color=[0,1,1];   % fimTotal - blue 
        CH(4).LineWidth = 3;
        CH(5).Color=[1,0,1];   % fimOpt - magenta
        CH(5).LineWidth = 3;
    end
end