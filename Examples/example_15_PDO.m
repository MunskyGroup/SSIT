%% SSIT/Examples/example_15_PDO

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
addpath(genpath('../'));
%addpath(genpath('tmpPropensityFunctions'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP
% example_6_SensitivityAnalysis
% example_7_FIM
% example_8_LoadingandFittingData_DataLoading
% example_9_LoadingandFittingData_MLE
% example_10_LoadingandFittingData_MHA

%% Load pre-run results:
% load('example_10_LoadingandFittingData_MHA.mat')

% View model summary:
STL1_4state_MH.summarizeModel

% Create a copy of the STL1 model for PDO:
STL1_4state_PDO = STL1_4state_MH;

%%
fimResults = STL1_4state_PDO.computeFIM(); 

% Get the number of cells using 'nCells':
cellCounts = STL1_4state_PDO.dataSet.nCells*ones(size(STL1_4state_PDO.tSpan));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Apply an Affine Poisson conditional probability distribution 
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

STL1_4state_PDO.plotMHResults(STL1_4state_MH_MHResults,[fimTotal],...
                              'log',[],figNew,plotColors)


% Find and store the total number of cells in your data set (already
% computed by SSIT when data was loaded in 
% example_8_LoadingandFittingData_DataLoading:
nTotal = sum(STL1_4state_PDO.dataSet.nCells);

%% Compute the optimal number of cells from the FIM results using the min. 
% inv determinant <x^{-1}> (all other parameters are known and fixed)
nCellsOpt = STL1_4state_PDO.optimizeCellCounts(fimResults,nTotal,...
                                                'Smallest Eigenvalue');
 
nCellsOptAvail = min(nCellsOpt,STL1_4state_PDO.dataSet.nCells)

fimOpt = STL1_4state_PDO.evaluateExperiment(fimResults,nCellsOpt,...
                                             diag(sig_log10.^2));

fimOptAvail = STL1_4state_PDO.evaluateExperiment(fimResults,...
                                        nCellsOptAvail,diag(sig_log10.^2));
figOpt = figure;
STL1_4state_PDO.plotMHResults(STL1_4state_MH_MHResults, [fimOpt,fimTotal],...
                              'log',[],figOpt,plotColors);
figOptAvail = figure;
STL1_4state_PDO.plotMHResults(STL1_4state_MH_MHResults,...
                   [fimOptAvail,fimTotal],'log',[],figOptAvail,plotColors);
 
f = figure;
set(f,'Position',[616   748   412   170])
bar([1:size(nCellsOpt,2)],nTotal,0.45)
hold on
bar([1:size(nCellsOpt,2)]+0.5,nCellsOpt,0.45)
set(gca,'xtick',[1:size(nCellsOpt,2)]+0.25,'xticklabel',...
    STL1_4state_PDO.dataSet.times,'fontsize',16,'ylim',[0,7000])
legend('Intuitive Design','Optimal Design')

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

STL1_4state_PDO_cyt = ...
 STL1_4state_PDO_cyt.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
        {'mRNA'}, {'RNA_STL1_total_TS3Full'}, {'RNA_STL1_cyto_TS3Full'},...
        'Binomial', true, [], {'Replica',1}, LegendLocation="northwest",...
        Title="4-state STL1 (PDO: Cytoplasmic mRNA)", FontSize=24,...
        XLabel="Total mRNA counts", YLabel="Cytoplasmic mRNA counts");

STL1_4state_PDO_nuc = ...
 STL1_4state_PDO_nuc.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
        {'mRNA'}, {'RNA_STL1_total_TS3Full'}, {'RNA_STL1_nuc_TS3Full'},...
        'Binomial', true, [], {'Replica',1}, LegendLocation="northwest",...
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
% Soon: Neural networks and hierarchical Bayes may be used to estimate λ
parGuess = [0, 1500, 5];

STL1_4state_PDO_intens = STL1_4state_PDO;
STL1_4state_PDO_intens = STL1_4state_PDO_intens.calibratePDO( ...
    'data/filtered_data_2M_NaCl_Step.csv', {'mRNA'},...
    {'RNA_STL1_total_TS3Full'}, {'STL1_avg_int_TS3Full'}, 'AffinePoiss',...
    true, parGuess, {'Replica',2}, LegendLocation="southeast", ...
    Title="4-state STL1 (PDO: Average intensity)", FontSize=24,...
    XLabel="True mRNA counts",YLabel="Average intensities (binned)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM + PDO analyses
%   * Ex(1): Analyze FIM with PDO for cytoplasmic & nuclear mRNA counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cytoplasm:
fimsPDO_cyt = STL1_4state_PDO_cyt.computeFIM([],'log');
fimPDO_cyt = STL1_4state_PDO_cyt.evaluateExperiment(fimsPDO_cyt,...
                                            nCellsOpt, diag(sig_log10.^2));

nCellsOptPDO_cyt = STL1_4state_PDO_cyt.optimizeCellCounts(...
                               fimsPDO_cyt, nTotal, 'Smallest Eigenvalue');
figPDO_cyt = figure;
STL1_4state_PDO_cyt.plotMHResults(STL1_4state_MH_MHResults,...
                         [fimPDO_cyt,fimTotal,fimOpt],'log',[],figPDO_cyt);

% Nucleus:
fimsPDO_nuc = STL1_4state_PDO_nuc.computeFIM([],'log');
fimPDO_nuc = STL1_4state_PDO_nuc.evaluateExperiment(fimsPDO_nuc,...
                                            nCellsOpt, diag(sig_log10.^2));

nCellsOptPDO_nuc = STL1_4state_PDO_nuc.optimizeCellCounts(...
                               fimsPDO_nuc, nTotal, 'Smallest Eigenvalue');
figPDO_nuc = figure;
STL1_4state_PDO_nuc.plotMHResults(STL1_4state_MH_MHResults,...
                         [fimPDO_nuc,fimTotal,fimOpt],'log',[],figPDO_nuc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   * Ex(2): Analyze FIM with PDO for average intensity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fimsPDOintens = STL1_4state_PDO_intens.computeFIM([],'log');
fimPDOintens = STL1_4state_PDO_intens.evaluateExperiment(fimsPDOintens,...
                                            nCellsOpt,diag(sig_log10.^2));

nCellsOptPDOintens = STL1_4state_PDO_intens.optimizeCellCounts(...
                               fimsPDOintens,nTotal,'Smallest Eigenvalue');

figintens = figure;
STL1_4state_PDO_intens.plotMHResults(STL1_4state_MH_MHResults,...
                        [fimPDOintens,fimTotal,fimOpt],'log',[],figintens);

%% Plot legend
axs = findall(figPDO_cyt, 'Type', 'axes');
ax  = axs(1);        
hold(ax,'on');

% Helpers:
near = @(c,tol,tgt) (numel(c)==3) && all(abs(c(:)'-tgt)<=tol);
isMagenta = @(c) near(c,0.15,[1 0 1]);
isCyan    = @(c) near(c,0.15,[0 1 1]);
isBlue    = @(c) near(c,0.15,[0 0 1]);
isGreen   = @(c) near(c,0.15,[0 1 0]);

% MCMC 90% credible interval (magenta dashed):
hMHell = findobj(ax,'Type','line','LineStyle','--');
hMHell = hMHell(arrayfun(@(h) isMagenta(h.Color), hMHell));

% FIM ellipses (solid lines):
hFIM = findobj(ax,'Type','line','LineStyle','-');

% Classify FIM ellipses by color:
hFIM_cyan  = hFIM(arrayfun(@(h) isCyan(h.Color),  hFIM));
hFIM_blue  = hFIM(arrayfun(@(h) isBlue(h.Color),  hFIM));
hFIM_green = hFIM(arrayfun(@(h) isGreen(h.Color), hFIM));

% Find MCMC samples (scatter) and MLE (square marker):
hSamples = findobj(ax,'Type','scatter');
if isempty(hSamples)
    cand = findobj(ax,'Type','line','Marker','o');
    hSamples = cand(~arrayfun(@(h) strcmp(get(h,'MarkerFaceColor'),'none'), cand));
end
hMLE = findobj(ax,'Type','line','Marker','s');

% Build legend in a sensible order:
L = []; names = {};
if ~isempty(hSamples),L(end+1)=hSamples(1);names{end+1}='MCMC samples';end
if ~isempty(hMLE),L(end+1)=hMLE(1); names{end+1}='MLE';end
if ~isempty(hMHell),L(end+1)=hMHell(1);names{end+1}='MCMC 90% CI';end
if ~isempty(hFIM_cyan),L(end+1)=hFIM_cyan(1);names{end+1}='FIM PDO';end
if ~isempty(hFIM_blue),L(end+1)=hFIM_blue(1);names{end+1}='FIM total';end
if ~isempty(hFIM_green),L(end+1)=hFIM_green(1);names{end+1}='FIM optimal';end

% Fallback: if color classification failed, just take first three FIM lines
if numel(L)<5
    remainingFIM = setdiff(hFIM, [hFIM_cyan; hFIM_blue; hFIM_green]);
    for k = 1:min(3, numel(remainingFIM))
        L(end+1) = remainingFIM(k);
        names{end+1} = sprintf('FIM #%d', k);
    end
end

lgd = legend(ax, L, names, 'Location','best');
lgd.FontSize = 12;


%% Save PDO models + results:
saveNames = unique({'STL1_4state_PDO'
    'fimTotal'
    'nTotal'
    'nCellsOpt'
    'nCellsOptAvail'
    'fimOpt'
    'fimOptAvail'
    'STL1_4state_PDO_cyt'
    'STL1_4state_PDO_nuc'
    'STL1_4state_PDO_intens'
    'fimsPDO_cyt'
    'fimsPDO_nuc'
    'fimsPDOintens'
    'fimPDO_cyt'
    'fimPDO_nuc'
    'fimPDOintens'
    'nCellsOptPDO_cyt'
    'nCellsOptPDO_nuc'
    'nCellsOptPDOintens'
    });
    
save('example_15_PDO',saveNames{:})
