%% SSIT/Examples/example_14_PDO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.4: Complex models
%   * Use probability distribution operators (PDOs) to handle distortion 
%     of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels, FSP solutions  
% from example_4_SolveSSITModels_FSP, sensitivities computed in 
% example_6_SensitivityAnalysis, FIM results from example_7_FIM,  
% loaded data from example_8_LoadingandFittingData_DataLoading, and
% Metropolis-Hastings results from example_10_LoadingandFittingData_MH
%clear
%close all

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_6_SensitivityAnalysis
% example_7_FIM
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE
% example_10_LoadingandFittingData_MH

%% Load pre-run results:
% load('example_7_FIM.mat')
% load('example_10_LoadingandFittingData_MH.mat')

% View model summary:
STL1_4state_MH.summarizeModel

% Create a copy of the STL1 model for PDO:
STL1_4state_PDO = STL1_4state_MH;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PDO Calculations
%% Ex(1): Calibrate PDO from cytoplasmic mRNA count data
% Calibrate Binomial PDOs from empirical data under the `missing spot'
% assumption from Vo et al. (2013).  'RNA_STL1_total_TS3Full' is the full 
% mRNA count and represents the `true' data; `RNA_STL1_cyto_TS3Full' is the
% mRNA count from only cytoplasm and `RNA_STL1_nuc_TS3Full' is the mRNA 
% count from only the nucleus.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STL1_4state_PDO_cyt = STL1_4state_PDO;
STL1_4state_PDO_nuc = STL1_4state_PDO;
parGuess = [];

STL1_4state_PDO_cyt = ...
 STL1_4state_PDO_cyt.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
      {'mRNA'}, {'RNA_STL1_total_TS3Full'}, {'RNA_STL1_cyto_TS3Full'},...
      'Binomial',true,parGuess,{'Replica',1},LegendLocation="northwest",...
       Title="4-state STL1 (PDO: Cytoplasmic mRNA)", FontSize=24,...
       XLabel="Total mRNA counts", YLabel="Cytoplasmic mRNA counts");

STL1_4state_PDO_nuc = ...
 STL1_4state_PDO_nuc.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
      {'mRNA'}, {'RNA_STL1_total_TS3Full'}, {'RNA_STL1_nuc_TS3Full'},...
      'Binomial',true,parGuess,{'Replica',1},LegendLocation="northwest",...
       Title="4-state STL1 (PDO: Nuclear mRNA)", FontSize=24,...
       XLabel="Total mRNA counts", YLabel="Nuclear mRNA counts");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Calibrate PDO from average intensity data
% Calibrate the PDO from empirical data. Here, the number of spots has
% been measured using different assays in data column 
% 'RNA_STL1_total_TS3Full' for the 'true' data set (mRNA spot counts) and 
% in the columns 'STL1_avg_int_TS3Full' for integrated intensity.  
% We calibrate an 'AffinePoiss' PDO where the obervation probability is a 
% Poisson distribution where the mean value is affine linearly related to 
% the true value: P(y|x) = Poiss(a0 + a1*x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make guesses for the PDO hyperparameters λ, in this case: λ₁, λ₂, λ₃
% Model: μ(x) = max(λ₁, λ₂ + λ₃·x)
parGuess = [0, 2000, 0];

STL1_4state_PDO_intens = STL1_4state_PDO;
STL1_4state_PDO_intens = STL1_4state_PDO_intens.calibratePDO( ...
    'data/filtered_data_2M_NaCl_Step.csv', {'mRNA'},...
    {'RNA_STL1_total_TS3Full'}, {'STL1_avg_int_TS3Full'}, 'AffinePoiss',...
    true, parGuess, {'Replica',1}, LegendLocation="southeast", ...
    Title="4-state STL1 (Affine PDO: Average intensity)", FontSize=24,...
    XLabel="True mRNA counts",YLabel="Average intensities (binned)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM + PDO analyses
%   * Analyze FIM with PDO for cytoplasmic & nuclear mRNA counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Free parameters:
freePars = 1:13;
sig_log10 = 0.1*ones(1,13);

% Get the number of cells using 'nCells' (same for both):
nTotal = sum(STL1_4state_PDO.dataSet.nCells);

% Cytoplasm:
fimsPDO_cyt = STL1_4state_PDO_cyt.computeFIM([],'log',[],freePars);
nCellsOptPDO_cyt = STL1_4state_PDO_cyt.optimizeCellCounts(fimsPDO_cyt,...
                                                nTotal, 'Trace');
fimPDO_cyt = STL1_4state_PDO_cyt.evaluateExperiment(fimsPDO_cyt,...
                                     nCellsOptPDO_cyt, diag(sig_log10.^2));

% Nucleus:
fimsPDO_nuc = STL1_4state_PDO_nuc.computeFIM([],'log',[],freePars);
nCellsOptPDO_nuc = STL1_4state_PDO_nuc.optimizeCellCounts(fimsPDO_nuc,...
                                                nTotal, 'Trace');
fimPDO_nuc = STL1_4state_PDO_nuc.evaluateExperiment(fimsPDO_nuc,...
                                     nCellsOptPDO_nuc, diag(sig_log10.^2));

% Plot the FIMs (cyt):
f20 = figure(20);
f21 = figure(21);
STL1_4state_PDO_cyt.plotFIMResults(fimPDO_cyt, 'log',...
    STL1_4state_PDO_cyt.parameters, PlotEllipses=true, EllipseFigure=f20,...
    EllipsePairs=[9 13; 4 13; 9 11; 10 11; 8 9; 12 13; 9 10; 10 12; 5 9],...
    FigureHandle=f21, Colors=struct('EllipseColors',[0.2 0.6 0.9],...
    'CenterSquare',[0.96,0.47,0.16]));

% Plot the FIMs (nuc):
f23 = figure(23);
STL1_4state_PDO_nuc.plotFIMResults(fimPDO_nuc, 'log',...
    STL1_4state_PDO_nuc.parameters(1:13), PlotEllipses=true,...
    EllipseFigure=f20, FigureHandle=f23,...
    EllipsePairs=[9 13; 4 13; 9 11; 10 11; 8 9; 12 13; 9 10; 10 12; 5 9],...
    Colors=struct('EllipseColors',[0 0 0],'CenterSquare',[0.96,0.47,0.16]));

% Plot the FIMs (free):
f24 = figure(24);
STL1_4state_FIM.plotFIMResults(STL1_4state_fimTotal_free, 'log',...
    STL1_4state_FIM.parameters(1:13), PlotEllipses=true, EllipseFigure=f20,...
    EllipsePairs=[9 13; 4 13; 9 11; 10 11; 8 9; 12 13; 9 10; 10 12; 5 9],...
    FigureHandle=f24, Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]));


%% Save PDO models + results:
saveNames = unique({'STL1_4state_PDO'
    'STL1_4state_PDO_cyt'
    'STL1_4state_PDO_nuc'
    'STL1_4state_PDO_intens'
    'fimsPDO_cyt'
    'fimsPDO_nuc'
    'fimPDO_cyt'
    'fimPDO_nuc'
    'nCellsOptPDO_cyt'
    'nCellsOptPDO_nuc'
    });
    
save('example_15_PDO',saveNames{:})
