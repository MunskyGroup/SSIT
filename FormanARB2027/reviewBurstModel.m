%% 
clear
clc
close all
addpath(genpath('..'))

%% Plots of means and fano factor versus parameters. Paper Figure 1

kr = 100;
gr = 1;
N = 100;
kArr = logspace(-2,3,N);
MEAN = zeros(N,N);
FANO = zeros(N,N);
for i = 1:N
    for j = 1:N
        F(i,j) = kArr(i)/(kArr(i)+kArr(j));
        MEAN(i,j) = kr/gr*F(i,j);
        FANO(i,j) = 1+((1-F(i,j))*kr)/(kArr(i)+kArr(j)+gr);
    end
end
figure(1); clf;
ax1 = axes;
[~,cF] = contourf(ax1,log10(kArr),log10(kArr),log10(FANO),linspace(0,log10(64),30)); hold on;
mu = 2;
g = mu*gr/kr;
kons = g*kArr/(1-g);
plot(log10(kArr), log10(kons), 'k--', 'LineWidth', 2)

mu = 25;
g = mu*gr/kr;
kons = g*kArr/(1-g);
plot(log10(kArr), log10(kons), 'k--', 'LineWidth', 2)

mu = 75;
g = mu*gr/kr;
kons = g*kArr/(1-g);
plot(log10(kArr), log10(kons), 'k--', 'LineWidth', 2)

colormap(ax1,jet)
cb = colorbar(ax1);
cb.Ticks = log10([1.0001,4,16,64]);
cb.TickLabels = {'1','4','16','64'};

set(gca,'xtick',[-2:3],'ytick',[-2:3],...
    'XTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'},...
    'YTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'},...
    'FontSize',16);
xlim([-2 3]);
ylim([-2 3]);

% ax2 = axes;
% [~,cM] = contour(ax2,log10(kArr),log10(kArr),MEAN,[2,25,75],'k--','LineWidth',2);
% ax2.Color = 'none';
% ax2.XColor = 'none';
% ax2.YColor = 'none';
% linkaxes([ax1 ax2]);
% ax2.Position = ax1.Position;



%% Plot key points on figure 1. Paper Figure 1
% find star
star_koff = 10;
mu = 2; 
g = (mu*gr)/kr;
star_kon = (g*star_koff)/(1-g);

% find collision point with mu == 25 line
% vary kon 
mu = 75;
g = (mu*gr)/kr;
final_kon = (g*star_koff)/(1-g);

% vary koff 
final_koff = (kr*star_kon)/(mu*gr)-star_kon;

% plot on figure 1
% Coordinates
x0 = log10(star_koff);
y0 = log10(star_kon);

x1 = log10(final_koff);
y1 = log10(star_kon);

x2 = log10(star_koff);
y2 = log10(final_kon);


% Squares
plot(ax1, x1, y1, 's', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b', ...
    'LineWidth', 3);

plot(ax1, x2, y2, 's', ...
    'MarkerSize', 12, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'r', ...
    'LineWidth', 3);

% Arrows
quiver(ax1, x0, y0, x1-x0, y1-y0, 0, ...
    Color='b', LineWidth=2, MaxHeadSize=0.25, ShowArrowHead='on');

quiver(ax1, x0, y0, x2-x0, y2-y0, 0, ...
    Color='r', LineWidth=2, MaxHeadSize=0.25, ShowArrowHead='on');

% Star
plot(ax1, x0, y0, 'k*', ...
    'MarkerSize', 15, 'LineWidth', 3);
% plot(ax1, x0, y0, 'k*', ...
%     'MarkerSize', 14, 'LineWidth', 1.5);

%% plot mean and std during transition from mu == 2 to mu == 75. Paper Figure 1
Model_chg = SSIT('Empty');
Model_chg.species = {'gON','mRNA'};
Model_chg.initialCondition = [0;0];
Model_chg.parameters = {'kon0', star_kon;...
    'koff0', star_koff;... 
    'kr',100;... 
    'g',gr;... 
    'kon1', star_kon; ...
    'koff1', star_koff; ...
    };

Model_chg.inputExpressions = {'I', ...
    't>=1'};

Model_chg = Model_chg.addReaction(struct(...
    'propensity',{'(kon0 + (kon1-kon0)*I)*(1-gON)'},...
    'stoichiometry',{{'gON',1}}));

Model_chg = Model_chg.addReaction(struct(...
    'propensity',{'(koff0 + (koff1-koff0)*I)*gON'},...
    'stoichiometry',{{'gON',-1}}));

Model_chg = Model_chg.addReaction(struct(...
    'propensity',{'kr*gON'},...
    'stoichiometry',{{'mRNA',1}}));

Model_chg = Model_chg.addReaction(struct(...
    'propensity',{'g*mRNA'},...
    'stoichiometry',{{'mRNA',-1}}));

Model_chg.fspOptions.initApproxSS = true;
Model_chg.tSpan = linspace(0,25,31);

Model_chg = Model_chg.solve;
Model_chg.plotFSP(plotType='marginals', SpeciesIdx=[2], indTimes=31, figureNums=2, Title='')

Model_chg.parameters{5,2} = star_kon;
Model_chg.parameters{6,2} = final_koff;
Model_chg = Model_chg.solve;
Model_chg.plotFSP(plotType='marginals', SpeciesIdx=[2], indTimes=31, figureNums=3, Title='')
Model_chg.plotFSP(plotType='meansAndDevs', SpeciesIdx=[2], Colors=[0 0 1], figureNums=5, Title='')

Model_chg.parameters{6,2} = star_koff;
Model_chg.parameters{5,2} = final_kon;
Model_chg = Model_chg.solve;
Model_chg.plotFSP(plotType='marginals', SpeciesIdx=[2], indTimes=31, figureNums=4, Title='')
Model_chg.plotFSP(plotType='meansAndDevs', SpeciesIdx=[2], Colors=[1 0 0], figureNums=5, Title='')

for figNum = [2 3 4]
    figure(figNum);
    xlim([0 120]);

    ax = gca;
    ax.Children.LineWidth = 4;
    ax.LineWidth = 2;
end

figure(5);
ax = gca;
set(findobj(ax, '-property', 'LineWidth'), 'LineWidth', 3);
ax.LineWidth = 2;

%% Export figures as SVG for PowerPoint
% Annual Review of Biochemistry
%
% Canvas: 6.33 x 7.9 inches
% Each figure row: 1.975 inches high

outputFolder = 'AnnualReview_Figures';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

fullWidth = 6.33;
thirdWidth = 6.33 / 3;
quarterHeight = 7.9 / 4;
sixtenthHeight = 7.9/16;

% Figure 1B
fig = figure(1);
ax = gca;

% Remove title and axis labels
ax.Title.String = '';
ax.XLabel.String = '';
ax.YLabel.String = '';
ax.XTickLabel = [];
ax.YTickLabel = [];
fig.Units = 'inches';
fig.Position(3:4) = [fullWidth 6*sixtenthHeight];

exportgraphics(fig, ...
    fullfile(outputFolder, 'figure1b.svg'), ...
    'ContentType', 'vector');

% Figure 1C-I, II, III
figNums = [2 3 4];
fileNames = {'figure1cI.svg', 'figure1cII.svg', 'figure1cIII.svg'};

for i = 1:3

    fig = figure(figNums(i));
    ax = gca;

    % Explicitly remove title
    sgtitle(fig, '');

    ax.Title.String = '';
    ax.Title.Visible = 'off';

    ax.Subtitle.String = '';
    ax.Subtitle.Visible = 'off';

    % Explicitly remove axis labels
    ax.XLabel.String = '';
    ax.XLabel.Visible = 'off';

    ax.YLabel.String = '';
    ax.YLabel.Visible = 'off';

    ax.YGrid = 'off';

    ax.XTickLabel = [];
    ax.YTickLabel = [];

    % Keep requested x-axis limits
    xlim(ax, [0 120]);

    % Set physical dimensions
    fig.Units = 'inches';
    fig.Position(3:4) = [thirdWidth 2*sixtenthHeight];

    % Export
    exportgraphics(fig, ...
        fullfile(outputFolder, fileNames{i}), ...
        'ContentType', 'vector');

end

% Figure 1D

fig = figure(5);
ax = gca;

% Remove title and axis labels
ax.Title.String = '';
ax.Title.Visible = 'off';

ax.XLabel.String = '';
ax.XLabel.Visible = 'off';

ax.YLabel.String = '';
ax.YLabel.Visible = 'off';

ax.XTickLabel = [];
ax.YTickLabel = [];

cb = colorbar(ax);
cb.TickLabels = [];

% Remove legend
leg = findobj(fig, 'Type', 'Legend');

if ~isempty(leg)
    delete(leg);
end

% Set physical dimensions
fig.Units = 'inches';
fig.Position(3:4) = [fullWidth quarterHeight];

% Export
exportgraphics(fig, ...
    fullfile(outputFolder, 'figure1d.svg'), ...
    'ContentType', 'vector');

disp('All SVG figures exported successfully.');



%% Simple experiment and eigenvector analysis
Model_chg.tSpan = linspace(0,7.5,31); % update time specific to the time scale
nCellsInExperiment = 0*Model_chg.tSpan;
nCellsInExperiment([1,31]) = 200;
FIMs = Model_chg.computeFIM(scale='log',freePars=[1:5],...
    observed={'mRNA'});
FIMTotal = Model_chg.totalFim(FIMs,nCellsInExperiment);

fprintf('FIM for step change in kon at 1 for intuitive design')
f = FIMTotal{1}
c = cond(f)
[V, D] = eig(f)
e = 1/2*(f*f')^-1
[V, D] = eig(e)

Model_chg.plotFIMResults(f, 'log', Model_chg.parameters(1:5), PlotEllipses=true, Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]))
% I really like this plot. It shows how the intial steady state has the
% least amount of information. 

%% Larger experiment and analysis
nCellsInExperiment = 0*Model_chg.tSpan;
nCellsInExperiment([1,3,5,10, 15, 20, 25, 31]) = 1000;
FIMs = Model_chg.computeFIM(scale='log',freePars=[1:5],...
    observed={'mRNA'});
FIMTotal = Model_chg.totalFim(FIMs,nCellsInExperiment);

fprintf('FIM for step change in kon at 1 for intuitive design')
f = FIMTotal{1}
c = cond(f) 
[V, D] = eig(f)
e = 1/2*(f*f')^-1
[V, D] = eig(e)

Model_chg.plotFIMResults(f, 'log', Model_chg.parameters(1:5),  PlotEllipses=true, Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]))
% see marginal improvement in intial steady state 


%% Setup
rng(1);

% True parameters
theta_true = [100, 50];

% Number of observations
n = 10;

% Generate data
ds = theta_true(1) + theta_true(2) * randn(n,1);

% Log-likelihood
logL = @(ds, m, s) sum(log(normpdf(ds, m, s)));

%% Single dataset

theta1_domain = linspace(0, theta_true(1)*2, 500);
likelihoods = zeros(size(theta1_domain));

for i = 1:length(theta1_domain)
    likelihoods(i) = logL(ds, theta1_domain(i), theta_true(2));
end

[value, idx] = max(likelihoods);
theta1_mle = theta1_domain(idx);

figure(10);
clf;

plot(theta1_domain, likelihoods, ...
    'Color', [0.15 0.35 0.65], ...
    'LineWidth', 2);

hold on;

% True parameter
xline(theta_true(1), ...
    '--k', ...
    'LineWidth', 2);

% MLE
xline(theta1_mle, ...
    '--r', ...
    'LineWidth', 2);

% Formatting
xlabel('\theta_1', 'FontSize', 14);
ylabel('Log-likelihood', 'FontSize', 14);

legend({'Log-likelihood', ...
        'True \theta_1', ...
        'MLE'}, ...
        'Location', 'best', ...
        'Box', 'off');

set(gca, ...
    'FontSize', 12, ...
    'LineWidth', 1.2, ...
    'TickDir', 'out', ...
    'Box', 'off');

grid off;

title(sprintf('Single dataset (n = %d)', n), ...
    'FontSize', 14, ...
    'FontWeight', 'normal');

%% Multiple datasets

all_thetas = [100, 50;
              100, 50;
              100, 100;
              100, 100];

N = [5, 10, 5, 10];

nDatasets = 1000;

for j = 1:length(N)

    theta_true = all_thetas(j,:);
    n = N(j);

    theta1_domain = linspace(0, theta_true(1)*2, 500);
    mle_estimates = zeros(nDatasets, 1);

    figure(10+j);
    clf;

    curveColor = [0.25 0.45 0.70];
    mleColor = [0.85 0.15 0.15];

    % Top: likelihood curves
    ax1 = axes('Position', [0.10 0.46 0.68 0.44]);
    hold(ax1, 'on');

    colors = ax1.ColorOrder;

    for k = 1:nDatasets

        % Generate dataset
        ds = theta_true(1) + theta_true(2) * randn(n,1);

        % Calculate log-likelihood
        likelihoods = zeros(size(theta1_domain));

        for i = 1:length(theta1_domain)
            likelihoods(i) = logL(ds, ...
                theta1_domain(i), ...
                theta_true(2));
        end

        % Exact MLE
        mle_estimates(k) = mean(ds);

        likelihoods = likelihoods - max(likelihoods);

        % Plot likelihood
        if k < 50
            c = colors(mod(k-1, size(colors,1))+1, :);

            plot(ax1, theta1_domain, likelihoods, ...
                'Color', c, ...
                'LineWidth', 1.5);
        end
    end

    ylim(ax1, [min(likelihoods(:)) 1]);

    % True parameter
    xline(ax1, theta_true(1), ...
        '--k', ...
        'LineWidth', 2);

    ax1.XAxisLocation = 'bottom';

    xlim(ax1, [0 theta_true(1)*2]);

    ylabel(ax1, 'Log-likelihood', ...
        'FontSize', 13);

    title(ax1, sprintf( ...
        '\\theta_1 = %g, \\theta_2 = %g, n = %d', ...
        theta_true(1), theta_true(2), n), ...
        'FontSize', 14, ...
        'FontWeight', 'normal');

    set(ax1, ...
        'FontSize', 11, ...
        'LineWidth', 1.2, ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'XColor', 'k');

    grid(ax1, 'off');

    % MLE markers
    plot(ax1, mle_estimates(1:50), ...
        zeros(size(mle_estimates(1:50))), ...
        'x', ...
        'Color', mleColor, ...
        'MarkerSize', 8, ...
        'LineWidth', 1.5);


    % Right: zoomed-in likelihood curves
    axZoom = axes('Position', [0.82 0.46 0.15 0.44]);
    hold(axZoom, 'on');

    colors = axZoom.ColorOrder;

    % Zoom window
    zoomWidth = 1;
    xlim(axZoom, ...
        [theta_true(1)-zoomWidth theta_true(1)+zoomWidth]);

    minLikelihoods = zeros(size(likelihoods));
    maxLikelihoods = zeros(size(likelihoods));
    zoom_theta_dom = linspace(theta_true(1)-zoomWidth, theta_true(1)+zoomWidth, length(theta1_domain));
    % Use same y limits as main panel
    ylim(axZoom, ax1.YLim);

    % Replot likelihood curves in zoomed region
    for k = 1:min(49, nDatasets)

        ds = theta_true(1) + theta_true(2) * randn(n,1);

        likelihoods = zeros(size(zoom_theta_dom));

        for i = 1:length(zoom_theta_dom)
            likelihoods(i) = logL(ds, ...
                zoom_theta_dom(i), ...
                theta_true(2));
        end

        likelihoods = likelihoods - max(likelihoods);
        minLikelihoods(k) = min(likelihoods);
        maxLikelihoods(k) = max(likelihoods);
        c = colors(mod(k-1, size(colors,1))+1, :);

        plot(axZoom, zoom_theta_dom, likelihoods, ...
            'Color', c, ...
            'LineWidth', 1.5);
    end

    ylim(axZoom, [-0.2 0]);

    % True parameter
    xline(axZoom, theta_true(1), ...
        '--k', ...
        'LineWidth', 2);

    % Formatting zoom panel
    set(axZoom, ...
        'FontSize', 10, ...
        'LineWidth', 1.2, ...
        'TickDir', 'out', ...
        'Box', 'on');

    grid(axZoom, 'off');

    xlabel(axZoom, '\theta_1');

    % Remove y-axis labels from zoom panel
    axZoom.YTickLabel = [];


    % Bottom: histogram
    ax2 = axes('Position', [0.10 0.08 0.68 0.38]);
    hold(ax2, 'on');

    histogram(ax2, mle_estimates, ...
        'NumBins', 20, ...
        'FaceColor', curveColor, ...
        'FaceAlpha', 0.75, ...
        'EdgeColor', 'none');

    % Flip histogram upside down
    ax2.YDir = 'reverse';

    % True parameter
    xline(ax2, theta_true(1), ...
        '--k', ...
        'LineWidth', 2);

    xlim(ax2, [0 theta_true(1)*2]);

    % Remove histogram x-axis
    ax2.XTick = [];
    ax2.XColor = 'none';

    % Remove histogram y-axis
    ax2.YTick = [];
    % ax2.YColor = 'none';

    ax2.Box = 'off';

    set(ax2, ...
        'FontSize', 11, ...
        'LineWidth', 1.2, ...
        'TickDir', 'out');

    grid(ax2, 'off');

    % Put likelihood axes in front
    uistack(ax1, 'top');

end












return 
%% Burst Frequency Model
Model = SSIT('Empty');

% compute kon0, s0, s1, and 

Model.species = {'gON','mRNA'};

Model.initialCondition = [0;0];

%                                     PRIOR
Model.parameters = {'kon0',0.01;...  % logn(-1,2)
    'koff0',0.1;...                   % logn(0,2)
    'kr',10;...                     % logn(1,2)
    'g',0.01;...                    % logn(-2,2)
    'kD',10;...                     % logn(1,2)
    'S0',1;...                      % NA (initial input concentration)
    'S1',5};                        % NA (final input concentration)

Model.inputExpressions = {'Iu','S0+(S1-S0)*(t>0)'};

Model = Model.addReaction(struct(...
    'propensity',{'kon0*(Iu/(kD+Iu))*(1-gON)'},...
    'stoichiometry',{{'gON',1}}));

Model = Model.addReaction(struct(...
    'propensity',{'koff0*gON'},...
    'stoichiometry',{{'gON',-1}}));

Model = Model.addReaction(struct(...
    'propensity',{'kr*gON'},...
    'stoichiometry',{{'mRNA',1}}));

Model = Model.addReaction(struct(...
    'propensity',{'g*mRNA'},...
    'stoichiometry',{{'mRNA',-1}}));

Model.fspOptions.initApproxSS = true;
Model.tSpan = linspace(0,300,31);

Model = Model.solve;

Model.plotFSP

%% 

%% Verification of FIM using CRLB (spread of MLE)
% nCellsInExperiment = 0*Model.tSpan;
% nCellsInExperiment([1,11,31]) = 200;
% nMLE = 200;
% Model.fittingOptions.modelVarsToFit = [1:5];
% MLE = Model.estimateMLEspread(nCells=nCellsInExperiment,observableSpecies={'mRNA'},nMLE=nMLE,simsSaveFile='BurstFIMSims.csv',freePars=[1:5],restart=true);
% MLE = Model.estimateMLEspread(nCells=nCellsInExperiment,observableSpecies={'mRNA'},nMLE=nMLE,simsSaveFile='BurstFIMSims.csv',freePars=[1:5],startPars=exp(MLE.mhSamples),restart=false);
% 
% FIMs = Model.computeFIM(scale='log',freePars=[1:5],...
%     observed={'mRNA'});
% FIMTotal = Model.totalFim(FIMs,nCellsInExperiment);
% 
% Model.plotMHResults(MLE,FIM=FIMTotal,fimScale='log',truncateChain=false);

%% 1D plots of likelihood function - Kr
% nCellsInExperiment = zeros(size(Model_chg.tSpan));
% nCellsInExperiment([2,13,31]) = 200;
% Model_chg.parameters{5,2} = star_kon;
% Model_chg.parameters{6,2} = final_koff;
% Model_chg.ssaOptions.Nexp = 1;
% Model_chg = Model_chg.solve;
% Model_chg.sampleDataFromFSP(saveFile='likelihoodData.csv',nCells=nCellsInExperiment,species2save={'mRNA'});
% Model_chg = Model_chg.loadData('likelihoodData.csv', {'mRNA', 'exp1_mRNA'});
% 
% pars = cell2mat(Model_chg.parameters(:,2));
% inter_idx = 35;
% par_idx = 5;
% 
% varyingPar = logspace(-2,2);
% likelihoods = zeros(size(varyingPar));
% 
% for i = 1:length(varyingPar)
%     pars(par_idx) = varyingPar(i); i
%     computeGrad = false;
%     if i == inter_idx
%         computeGrad = true;
%     end
%     [logL, grad, ~] = Model_chg.computeLikelihood(pars, [], computeGrad);
%     likelihoods(i) = logL;
%     if i == inter_idx
%         gradients = grad;
%     end
% end
% 
% fims = Model_chg.computeFIM(scale='log',freePars=[1:6],...
%     observed={'mRNA'});
% 
% nCellsInExperiment = Model.dataSet.nCells;
% fim = Model_chg.totalFim(fims,nCellsInExperiment);
% 
% %%
% figure(6)
% semilogx(varyingPar, likelihoods, 'b-'); hold on
% grid on
% 
% % Save original axis limits
% xlim0 = xlim;
% ylim0 = ylim;
% 
% % Point and slope
% x1 = varyingPar(inter_idx);
% y1 = likelihoods(inter_idx);
% m = gradients(par_idx)*x1*log(10);
% 
% % ---------------- Tangent line ----------------
% lineDecade = 0.5;
% x = logspace(log10(x1)-lineDecade, log10(x1)+lineDecade, 50);
% y = y1 + m*(log10(x)-log10(x1));
% semilogx(x, y, 'Color', [0 0.5 0],  'LineWidth', 2)
% 
% % Point
% semilogx(x1, y1, 'ro', 'MarkerFaceColor', 'r')
% 
% % % ---------------- Small slope triangle ----------------
% % triDecade = 0.08;
% % 
% % x2 = x1*10^triDecade;
% % y2 = y1;
% % y3 = y1 + m*triDecade;
% % 
% % % Run
% % semilogx([x1 x2], [y1 y2], 'k-', 'LineWidth', 1.5)
% % 
% % % Rise
% % semilogx([x2 x2], [y2 y3], 'k-', 'LineWidth', 1.5)
% 
% % % ---------------- Restore original axes ----------------
% % xlim(xlim0)
% % ylim(ylim0)
% 
% % Peak parameter and likelihood
% pars = cell2mat(Model_chg.parameters(:,2));
% [~, idx] = min(abs(varyingPar - pars(par_idx)));
% x1 = varyingPar(idx);
% y1 = likelihoods(idx);
% 
% % Hessian / curvature
% d2 = -fim{1}(par_idx,par_idx);
% 
% % Work entirely in log10(parameter) space
% u1 = log10(x1);
% 
% % Small region around peak
% width = 0.25;
% u = linspace(u1-width, u1+width, 100);
% x_quad = 10.^u;
% 
% % Convert curvature from linear parameter space to log space
% d2_log = d2 * (x1)^2;
% 
% % Taylor expansion around peak (first derivative ~ 0)
% y_quad = y1 + 0.5*d2_log*(u-u1).^2;
% 
% % Plot quadratic
% semilogx(x_quad, y_quad, ...
%     'Color',[1 0.5 0], ...
%     'LineWidth',2);
% 
% xlabel('varyingPar')
% ylabel('Likelihood')
% title(sprintf('Slope = %.3g', m))
% legend('Likelihood', 'Tangent', 'Point', 'Location', 'best')
% 
% 
% return
% 
% %% 1D plots of likelihood function - gamma
% pars = cell2mat(Model.parameters(:,2));
% varyingPar = logspace(-4,0);
% likelihoods = zeros(size(varyingPar));
% gradients = cell(size(varyingPar));
% 
% for i = 1:length(varyingPar)
%     pars(4) = varyingPar(i);
%     computeGrad = false;
%     if i == inter_idx
%         computeGrad = true;
%     end
%     [logL, grad, ~] = Model.computeLikelihood(pars, [], true);
%     likelihoods(i) = logL;
%     if i == inter_idx
%         gradients = grad;
%     end
% end
% 
% figure(4)
% semilogx(varyingPar, likelihoods)


%% FIM Calculations
Sarray = [1:5];
Model.tSpan = [0:30];
Model.solutionScheme = 'fspsens';
FIM = cell(length(Sarray),length(Sarray),length(Model.tSpan));
for iS0 = 1:length(Sarray)
    for iS1 = 1:length(Sarray)
        Model = Model.changeParameter({'S0',Sarray(iS0);'S1',Sarray(iS1)-Sarray(iS0)});
        Model = Model.solve;
        FIM(iS0,iS1,:) = Model.computeFIM(freePars=(1:4),scale='log');
    end
end

%% FIM for different experiment designs.
Ncells = 600;
% Measurment at one steady state values.
iS = 3;
FIM_One_SS = Ncells*FIM{iS,iS,1};
disp(['Determinant of FIM for one SS measurement: ',num2str(det(FIM_One_SS))])

% Measurment at two steady state values.
iS1 = 1;
iS2 = 5;
FIM_Two_SS = Ncells/2*FIM{iS1,iS1,1}+Ncells/2*FIM{iS2,iS2,1};
disp(['Determinant of FIM for two SS measurements: ',num2str(det(FIM_Two_SS))])

% Measurement at change from one SS to another at three time points.
iS1 = 1;
iS2 = 5;
itimes = [1,11,31];
FIM_Dynamic = 0;
for it = 1:length(itimes)
    FIM_Dynamic = FIM_Dynamic + Ncells/length(itimes)*FIM{iS1,iS2,itimes(it)};
end
disp(['Determinant of FIM for dynamic measurements: ',num2str(det(FIM_Dynamic))])

%% Plots of FIM predicted uncertainties
freePars = [1:4];
f1 = figure(11);
f2 = figure(12);
f3 = figure(13);

% The following plots the heatmap showing the
Model.plotFIMResults(FIM_Dynamic^(-1)/log(10)^2, 'log',...
    Model.parameters(freePars,1),...
    [Model.parameters{freePars,2}],...
    PlotEllipses=true,EllipseFigure=f1,...
    FigureHandle=f3,...
    Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]),...
    LogThreshold=-1,...
    HeatMapType='invfim',...
    MatrixType='invfim');


%% Optimized Experiment Design
% AllFims = reshape(FIM,numel(FIM),1);
allFims = {};%cell(numel(FIM),1);
indsFims = [];zeros(numel(FIM),3);
k = 0;
for iS0 = 1:length(Sarray)
    for iS1 = 1:length(Sarray)
        for iT = 1:length(Model.tSpan)
            k = k+1;
            allFims(k,1) = FIM(iS0,iS1,iT);
            indsFims(k,:) = [iS0,iS1,iT];
        end
    end
end
OptExperiment = Model.optimizeCellCounts(allFims,600,'D-opt');
J = find(OptExperiment);
disp(['Optimized Experiment Design:'])
for j = 1:length(J)
    % Store optimized parameters and their corresponding indices
    optimizedParams = OptExperiment(J(j));
    paramIndices = indsFims(J(j), :);
    if paramIndices(3)==1||paramIndices(1)==paramIndices(2) % SS experiment
        disp(['   ',num2str(optimizedParams),' cells at steady state for S0 = ',num2str(Sarray(paramIndices(1)))])
    else
        disp(['   ',num2str(optimizedParams),' cells at time ',num2str(Model.tSpan(paramIndices(3))),' for S0 = ',num2str(Sarray(paramIndices(1))),' and S1 = ',num2str(Sarray(paramIndices(2)))])
    end
end

FIM_Opt = 0;
for i = 1:length(OptExperiment)
    FIM_Opt = FIM_Opt + OptExperiment(i)*allFims{i};
end
disp(['Determinant of FIM for optimized measurements: ',num2str(det(FIM_Opt))])

return
%% Generalized model
% 
% In the most general form of this model, the parameter 'alpha' determines
% if the model is burst size regulate (koff, alpha=0) or burst frequency
% regulated (kon, alpha=1) or something in between.
ModelGen = SSIT('Empty');

ModelGen.species = {'gON','mRNA'};

ModelGen.initialCondition = [0;0];

%                                     PRIOR
ModelGen.parameters = {'kon0',0.1;...  % logn(-1,2)
    'koff0',1;...                   % logn(0,2)
    'kr',10;...                     % logn(1,2)
    'g',0.01;...                    % logn(-2,2)
    'kD',10;...                     % logn(1,2)
    'alpha',1e-6;...                   % n(0,2)
    'S0',1;...                      % NA (initial input concentration)
    'S1',5};                        % NA (final input concentration)

ModelGen.inputExpressions = {'Iu','S0+S1*(t>0)'};

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'kon0*((1-alpha)+(alpha)*Iu/(kD+Iu))*(1-gON)'},...
    'stoichiometry',{{'gON',1}}));

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'koff0*((alpha)+(1-alpha)*(kD+Iu)/Iu)*gON'},...
    'stoichiometry',{{'gON',-1}}));

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'kr*gON'},...
    'stoichiometry',{{'mRNA',1}}));

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'g*mRNA'},...
    'stoichiometry',{{'mRNA',-1}}));

ModelGen.tSpan = linspace(0,300,31);

ModelGen = ModelGen.solve;

ModelGen.plotFSP

%% Parameters for key points.
kr = 100;
gr = 1;
KD = 10;
kon0a = 2.5;
koff0a = 10;
S0 = KD/24;
S1 = 23*KD/24;

%                                     PRIOR
ModelGen.parameters = {'kon0',kon0;...  % logn(-1,2)
    'koff0',koff0;...                   % logn(0,2)
    'kr',kr;...                     % logn(1,2)
    'g',gr;...                    % logn(-2,2)
    'kD',KD;...                     % logn(1,2)
    'alpha',-6;...                   % n(0,2)
    'S0',S0;...                      % NA (initial input concentration)
    'S1',S1};                        % NA (final input concentration)

ModelGen.tSpan = linspace(0,6,61);

ModelGen.fspOptions.initApproxSS = true;
ModelGen = ModelGen.solve;

ModelGen.plotFSP
