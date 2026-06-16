function [results, figHandle, ax] = plotSeqModelSummary(dataFileName, modelDir, varargin)
% plotSeqModelSummary
% Collects fitted SSIT gene model results from .mat files and makes a 2x2
% summary figure similar to the end-of-script plotting block you shared.
%
% USAGE
%   [results, figH] = plotSeqModelSummary( ...
%       'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv', ...
%       'seqModels_fit', ...
%       'Title', 'Dex response summary', ...
%       'FigureNumber', 1);
%
% INPUTS
%   dataFileName : (char/string) CSV table used only to discover gene names.
%   modelDir     : (char/string) folder containing "Model_<GENE>.mat" files.
%
% NAME-VALUE OPTIONS (all optional)
%   'GeneColumns'          : indices of columns in the CSV that correspond to genes
%                            (default mimics your logic: 2 : (end-4))
%   'TimeScale'            : scalar multiplier used in your KON/KOFF transforms (default 0.956)
%   'FigureNumber'         : figure number to use (default 1)
%   'Title'                : overall figure title (default "")
%   'KonLabel'             : x-label for KON in subplots 2-4 (default 'K_{ON}')
%   'KoffLabel'            : y-label for KOFF in subplots 2-4 (default 'K_{OFF}')
%   'RatioKonLabel'        : x-label for ratio plot (default 'max(K_{ON}) / min(K_{ON})')
%   'RatioKoffLabel'       : y-label for ratio plot (default 'max(K_{OFF}) / min(K_{OFF})')
%   'ColorbarLabel'        : colorbar label (default 'Mean Fold Change')
%   'FontSize'             : axes font size (default 16)
%   'KonKoffXLim'          : [min max] for KON in subplots 2-4 (default [1e-6 1e2])
%   'KonKoffYLim'          : [min max] for KOFF in subplots 2-4 (default [1e-4 1e4])
%   'KonKoffXTicks'        : ticks for KON in subplots 2-4 (default 10.^(-6:4:2))
%   'KonKoffYTicks'        : ticks for KOFF in subplots 2-4 (default 10.^(-4:4:4))
%   'RatioXLim'            : [min max] for ratio plot xlim (default [1 1e6])
%   'RatioXTicks'          : ticks for ratio plot x-axis (default 10.^(0:2:6))
%   'ColorLimits'          : [cmin cmax] for colorbar (default [0 9])
%   'ColorTicks'           : ticks for colorbar (default 0:2:8)
%   'ColorTickLabels'      : cellstr labels for color ticks (default {'2^0','2^2',...})
%
% OUTPUTS
%   results   : struct with geneNames, parameters, summary metrics, and group labels J.
%   figHandle : handle to the figure
%   ax        : struct of axes handles (ax.ratio, ax.group2, ax.group3, ax.group4)

% ----------------------------
% Parse inputs / options
% ----------------------------
p = inputParser;
p.addRequired('dataFileName', @(s) ischar(s) || isstring(s));
p.addRequired('modelDir',     @(s) ischar(s) || isstring(s));

p.addParameter('GeneColumns', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('TimeScale', 0.956, @(x) isnumeric(x) && isscalar(x) && isfinite(x));

p.addParameter('FigureNumber', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Title', "", @(s) ischar(s) || isstring(s));

p.addParameter('KonLabel', 'K_{ON}', @(s) ischar(s) || isstring(s));
p.addParameter('KoffLabel','K_{OFF}',@(s) ischar(s) || isstring(s));
p.addParameter('RatioKonLabel','max(K_{ON}) / min(K_{ON})', @(s) ischar(s) || isstring(s));
p.addParameter('RatioKoffLabel','max(K_{OFF}) / min(K_{OFF})', @(s) ischar(s) || isstring(s));
p.addParameter('ColorbarLabel','Mean Fold Change', @(s) ischar(s) || isstring(s));

p.addParameter('FontSize', 16, @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('KonKoffXLim', [1e-6 1e2], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('KonKoffYLim', [1e-4 1e4], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('KonKoffXTicks', 10.^(-6:4:2), @(x) isnumeric(x) && isvector(x));
p.addParameter('KonKoffYTicks', 10.^(-4:4:4), @(x) isnumeric(x) && isvector(x));

p.addParameter('RatioXLim', [1 1e6], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('RatioXTicks', 10.^(0:2:6), @(x) isnumeric(x) && isvector(x));

p.addParameter('ColorLimits', [0 9], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('ColorTicks', 0:2:8, @(x) isnumeric(x) && isvector(x));
p.addParameter('ColorTickLabels', {'2^0','2^2','2^4','2^6','2^8'}, @(c) iscellstr(c) || isstring(c));

p.parse(dataFileName, modelDir, varargin{:});
opt = p.Results;

% ----------------------------
% Read table and determine gene list
% ----------------------------
TAB = readtable(opt.dataFileName);

% Variable names are the table column headers (what you actually want here).
varNames = TAB.Properties.VariableNames;

% Default gene column selection: mimic your "2:end-4" behavior.
if isempty(opt.GeneColumns)
    geneCols = 2:(numel(varNames)-4);
else
    geneCols = opt.GeneColumns;
end

geneNames = varNames(geneCols);
nGenes = numel(geneNames);

% ----------------------------
% Preallocate arrays for speed + clarity
% ----------------------------
% NOTE: You used 9 parameters; keep as-is.
Pars        = nan(nGenes, 9);

Mean1       = nan(1, nGenes);
Var1        = nan(1, nGenes);
Fano1       = nan(1, nGenes);
CV1         = nan(1, nGenes);

MeanEnd     = nan(1, nGenes);
VarEnd      = nan(1, nGenes);
FanoEnd     = nan(1, nGenes);
CVEnd       = nan(1, nGenes);

DeltaFano   = nan(1, nGenes);
MaxDeltaFano= nan(1, nGenes);   % (bug fix) was used in your snippet but not preallocated
DeltaMean   = nan(1, nGenes);
MaxDeltaMean= nan(1, nGenes);
DeltaVar    = nan(1, nGenes);
DeltaCV     = nan(1, nGenes);

fileExists  = false(1, nGenes);

% Small epsilon to avoid division-by-zero or log(0) issues.
epsMean = 1e-6;

% ----------------------------
% Load each model file and extract summary metrics
% ----------------------------
for iGene = 1:nGenes
    modelName = ['Model_' geneNames{iGene}];
    matPath   = fullfile(opt.modelDir, [modelName '.mat']);

    if ~exist(matPath, 'file')
        % Leave NaNs in place for missing genes; keep a boolean flag.
        fileExists(iGene) = false;
        continue
    end

    fileExists(iGene) = true;

    % Load *without* eval: the .mat contains a variable named like "Model_GENE".
    S = load(matPath);

    % Guard: if the expected variable isn't present, mark missing.
    if ~isfield(S, modelName)
        fileExists(iGene) = false;
        continue
    end

    m = S.(modelName);

    % ---- Parameters ----
    % m.parameters is usually an Nx2 cell array: {name, value}.
    % Convert the value column into a numeric vector.
    try
        parVals = cellfun(@(v) double(v), m.parameters(:,2));
    catch
        % Fallback if values aren't directly double-castable
        parVals = cell2mat(m.parameters(:,2));
    end

    % Store first 9 (or pad/truncate safely).
    Pars(iGene, 1:min(9,numel(parVals))) = parVals(1:min(9,numel(parVals)));

    % ---- Dataset stats ----
    % Assumes m.dataSet.mean and m.dataSet.var exist and are vectors over time.
    mu  = m.dataSet.mean(:)';
    va  = m.dataSet.var(:)';

    % First time point
    Mean1(iGene) = mu(1) + epsMean;
    Var1(iGene)  = va(1);
    Fano1(iGene) = Var1(iGene) / Mean1(iGene);
    CV1(iGene)   = Var1(iGene) / (Mean1(iGene)^2);

    % Last time point
    MeanEnd(iGene) = mu(end) + epsMean;
    VarEnd(iGene)  = va(end);
    FanoEnd(iGene) = VarEnd(iGene) / MeanEnd(iGene);
    CVEnd(iGene)   = VarEnd(iGene) / (MeanEnd(iGene)^2);

    % Fold / delta summaries (relative to first time point)
    DeltaFano(iGene)    = (va(end)/ (mu(end)+epsMean)) - (va(1)/ (mu(1)+epsMean));
    MaxDeltaFano(iGene) = max(va ./ (mu+epsMean))      - (va(1)/ (mu(1)+epsMean));

    DeltaCV(iGene)      = (va(end)/ (mu(end)+epsMean)^2) - (va(1)/ (mu(1)+epsMean)^2);

    DeltaMean(iGene)    = (mu(end)+epsMean) / (mu(1)+epsMean);
    MaxDeltaMean(iGene) = max(mu+epsMean)    / (mu(1)+epsMean);

    DeltaVar(iGene)     = (va(end)+epsMean) / (va(1)+epsMean);
end

% ----------------------------
% Compute KON/KOFF transforms you used for plotting
% ----------------------------
% Your code used Pars(:,3), Pars(:,4), Pars(:,5), Pars(:,6).
% Keep the same formula; interpret opt.TimeScale as the 0.956 constant.
ts = opt.TimeScale;

Xvals  = (Pars(:,3) + Pars(:,4)*ts) ./ Pars(:,3);  % ratio on KON-related terms
Yvals  = (1 + Pars(:,6)*ts);                       % ratio on KOFF-related terms

Xvals1 = Pars(:,3);
Xvals2 = (Pars(:,3) + Pars(:,4)*ts);

Yvals1 = Pars(:,5);
Yvals2 = Pars(:,5) ./ (1 + Pars(:,6)*ts);

% ----------------------------
% Make figure
% ----------------------------
figHandle = figure(opt.FigureNumber); clf(figHandle);

% Axes handle struct (so caller can tweak later if desired)
ax = struct();

% J labels which group each gene into one of three categories
J = zeros(size(Xvals1));  % 1,2,3 correspond to subplots 3,2,4 in your logic

% --- Subplots 2-4: line segments from (min)->(max) for each gene ---
for i = 1:nGenes
    % Skip genes without files or with missing parameters
    if ~fileExists(i) || any(isnan([Xvals1(i),Xvals2(i),Yvals1(i),Yvals2(i)]))
        continue
    end

    % Your decision rule: compare KON ratio vs KOFF ratio thresholds
    konRatio  = Xvals2(i)/Xvals1(i);
    koffRatio = Yvals1(i)/Yvals2(i);

    if konRatio >= 3*koffRatio
        ax.group3 = subplot(2,2,3);
        plot([Xvals1(i), Xvals2(i)], [Yvals1(i), Yvals2(i)], '-s'); hold on
        J(i) = 1;
    elseif konRatio <= (1/3)*koffRatio
        ax.group2 = subplot(2,2,2);
        plot([Xvals1(i), Xvals2(i)], [Yvals1(i), Yvals2(i)], '-o'); hold on
        J(i) = 2;
    else
        ax.group4 = subplot(2,2,4);
        plot([Xvals1(i), Xvals2(i)], [Yvals1(i), Yvals2(i)], '-^'); hold on
        J(i) = 3;
    end

    % Apply consistent formatting to whichever subplot we're currently in
    set(gca, ...
        'XScale','log','YScale','log', ...
        'FontSize', opt.FontSize, ...
        'YLim', opt.KonKoffYLim, ...
        'XLim', opt.KonKoffXLim, ...
        'XTick', opt.KonKoffXTicks, ...
        'YTick', opt.KonKoffYTicks);

    xlabel(opt.KonLabel);
    ylabel(opt.KoffLabel);
end

% --- Subplot 1: ratio scatter colored by log2(MaxDeltaMean) ---
ax.ratio = subplot(2,2,1); hold off

% Color quantity: log2(MaxDeltaMean) + eps to avoid log2(0)
cVal = log2(MaxDeltaMean + epsMean);

% Plot each category with different markers (as you did)
scatter(Xvals(J==1), Yvals(J==1), 100, cVal(J==1), 's', 'filled'); hold on
scatter(Xvals(J==2), Yvals(J==2), 100, cVal(J==2), 'o', 'filled'); hold on
scatter(Xvals(J==3), Yvals(J==3), 100, cVal(J==3), '^', 'filled'); hold on

set(gca, ...
    'XScale','log','YScale','log', ...
    'FontSize', opt.FontSize, ...
    'XTick', opt.RatioXTicks, ...
    'XLim',  opt.RatioXLim);

% Reference dashed lines at ratio=2 (as in your code)
xL = get(gca,'XLim');
yL = get(gca,'YLim');
plot(xL, [2 2], 'k--');
plot([2 2], yL, 'k--');

% Colorbar formatting
C = colorbar;
C.Limits     = opt.ColorLimits;
C.Ticks      = opt.ColorTicks;
C.TickLabels = cellstr(opt.ColorTickLabels);
C.Label.String = opt.ColorbarLabel;

xlabel(opt.RatioKonLabel);
ylabel(opt.RatioKoffLabel);

% Optional overall title
if strlength(string(opt.Title)) > 0
    sgtitle(string(opt.Title), 'FontSize', opt.FontSize + 2);
end

% ----------------------------
% Package results for caller
% ----------------------------
results = struct();
results.dataFileName  = string(opt.dataFileName);
results.modelDir      = string(opt.modelDir);
results.geneNames     = geneNames;
results.fileExists    = fileExists;

results.Pars          = Pars;

results.Mean1         = Mean1;
results.Var1          = Var1;
results.Fano1         = Fano1;
results.CV1           = CV1;

results.MeanEnd       = MeanEnd;
results.VarEnd        = VarEnd;
results.FanoEnd       = FanoEnd;
results.CVEnd         = CVEnd;

results.DeltaFano     = DeltaFano;
results.MaxDeltaFano  = MaxDeltaFano;
results.DeltaMean     = DeltaMean;
results.MaxDeltaMean  = MaxDeltaMean;
results.DeltaVar      = DeltaVar;
results.DeltaCV       = DeltaCV;

results.Xvals         = Xvals;
results.Yvals         = Yvals;
results.J             = J;  % grouping label used for markers/categories

end