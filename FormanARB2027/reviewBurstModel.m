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

%% Export figures for Paper Figure 1
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
nCellsInExperiment([1]) = 1;
FIMs = Model_chg.computeFIM(scale='log',freePars=[1:4],...
    observed={'mRNA'});
FIMTotal = Model_chg.totalFim(FIMs,nCellsInExperiment);

fprintf('FIM for step change in kon at 1 for intuitive design')
f = FIMTotal{1}
c = cond(f)
[V, D] = eig(f)
e = 1/2*(f*f')^-1
[V, D] = eig(e)
[lambda,idx] = sort(diag(D),'descend');
D = diag(lambda);
V = V(:,idx);

% Model_chg.plotFIMResults(f, 'log', Model_chg.parameters(1:4), PlotEllipses=true, Colors=struct('EllipseColors',[0.9 0.6 0.2],...
%     'CenterSquare',[0.96,0.47,0.16]))

figure(18)
plotHeatmap( ...
    f, ...
    Model_chg.parameters(1:4,1), ...
    Model_chg.parameters(1:4,1), ...
    'FIM');


figure(19)
plotHeatmap( ...
    e, ...
    Model_chg.parameters(1:4,1), ...
    Model_chg.parameters(1:4,1), ...
    'FIM inverse');

figure(20)
plotHeatmap( ...
    V * D, ...
    Model_chg.parameters(1:4,1), ...
    {'EV1', 'EV2', 'EV3', 'EV4'}, ...
    'Eigen Vectors of FIM inverse');


% I really like this plot. It shows how the intial steady state has the
% least amount of information. 

%% Larger experiment and analysis
% nCellsInExperiment = 0*Model_chg.tSpan;
% nCellsInExperiment([1,3,5,10, 15, 20, 25, 31]) = 1000;
% FIMs = Model_chg.computeFIM(scale='log',freePars=[1:5],...
%     observed={'mRNA'});
% FIMTotal = Model_chg.totalFim(FIMs,nCellsInExperiment);

% fprintf('FIM for step change in kon at 1 for intuitive design')
% f = FIMTotal{1}
% c = cond(f) 
% [V, D] = eig(f)
% e = 1/2*(f*f')^-1
% [V, D] = eig(e)
% 
% Model_chg.plotFIMResults(f, 'log', Model_chg.parameters(1:5), [Model_chg.parameters{1:5,2}] ,PlotEllipses=true, Colors=struct('EllipseColors',[0.9 0.6 0.2],...
%     'CenterSquare',[0.96,0.47,0.16]))
% see marginal improvement in intial steady state 


%% Setup - MLE FIM relationship - Gaussian
rng(1);

% True parameters
theta_true = [100, 50];

% Number of observations
n = 10;

% Generate data
ds = theta_true(1) + theta_true(2) * randn(n,1);

% Log-likelihood
logL = @(ds, m, s) sum(log(normpdf(ds, m, s)));

%% MLE FIM relationship - Single dataset - Gaussian

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

%% MLE FIM relationship - Multiple datasets - Gaussian

all_thetas = [100, 50;
              100, 50;
              100, 100;
              100, 100];

N = [10, 5, 10, 5];

nDatasets = 10000;

% for j = 1:length(N)
for j = 1:1

    theta_true = all_thetas(j,:);
    n = N(j);

    theta1_domain = linspace(0, theta_true(1)*2, 1000);
    mle_estimates = zeros(nDatasets, 1);

    % Store likelihood curves so the zoom is exactly the same data
    all_likelihoods = zeros(nDatasets, length(theta1_domain));

    figure(10+j);
    clf;

    curveColor = [0.25 0.45 0.70];
    mleColor = [0.85 0.15 0.15];

    % Use one color order for both panels
    colors = turbo(49);

    % Top: likelihood curves
    ax1 = axes('Position', [0.10 0.46 0.68 0.44]);
    hold(ax1, 'on');

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

        % Shift maximum to zero
        likelihoods = likelihoods - max(likelihoods);

        % Store likelihood
        all_likelihoods(k,:) = likelihoods;

        % Plot first 49 datasets
        if k < 50
            plot(ax1, theta1_domain, likelihoods, ...
                'Color', colors(k,:), ...
                'LineWidth', 1);
        end
    end

    ylim(ax1, [min(all_likelihoods(1:49,:), [], 'all') 0.5]);

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

    % MLE markers directly on x-axis
    plot(ax1, mle_estimates(1:50), ...
        ones(size(mle_estimates(1:50)))*min(all_likelihoods(1:49,:), [], 'all'), ...
        'x', ...
        'Color', mleColor, ...
        'MarkerSize', 8, ...
        'LineWidth', 1);


    % Right: zoomed-in likelihood curves
    % zoomWidth = 1;
    % 
    % leftIdx = find(theta1_domain >= theta_true(1)-zoomWidth, 1);
    % rightIdx = find(theta1_domain <= theta_true(1)+zoomWidth, 1, 'last');
    % 
    % axZoom = axes('Position', [0.82 0.46 0.15 0.44]);
    % hold(axZoom, 'on');
    % 
    % xlim(axZoom, ...
    %     [theta_true(1)-zoomWidth theta_true(1)+zoomWidth]);
    % 
    % % Use same y range as main plot
    % % ylim(axZoom, [min(all_likelihoods(1:49, leftIdx:rightIdx), [], 'all'), 0.01]);
    % ylim(axZoom, [-3, 0.01]);
    % % Plot EXACT SAME curves with EXACT SAME colors
    % for k = 1:49
    %     plot(axZoom, theta1_domain, all_likelihoods(k,:), ...
    %         'Color', colors(k,:), ...
    %         'LineWidth', 1.5);
    % end
    % 
    % % True parameter
    % xline(axZoom, theta_true(1), ...
    %     '--k', ...
    %     'LineWidth', 2);
    % 
    % % Formatting zoom panel
    % set(axZoom, ...
    %     'FontSize', 10, ...
    %     'LineWidth', 1.2, ...
    %     'TickDir', 'out', ...
    %     'Box', 'on');
    % 
    % grid(axZoom, 'off');
    % 
    % xlabel(axZoom, '\theta_1');
    % 
    % % Remove y-axis labels
    % axZoom.YTickLabel = [];


    % Bottom: histogram
    ax2 = axes('Position', [0.10 0.08 0.68 0.38]);
    hold(ax2, 'on');

    histogram(ax2, mle_estimates, ...
        'NumBins', 50, ...
        'FaceColor', curveColor, ...
        'FaceAlpha', 0.75, ...
        'EdgeColor', 'none', ...
        'Normalization', 'pdf');

    pmle = normpdf(theta1_domain, theta_true(1), theta_true(2)/sqrt(n));
    plot(ax2, theta1_domain, pmle)

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

    % Remove histogram y-axis numbers
    ax2.YTick = [];

    ax2.Box = 'off';

    set(ax2, ...
        'FontSize', 11, ...
        'LineWidth', 1.2, ...
        'TickDir', 'out');

    grid(ax2, 'off');

    % Put likelihood axes in front
    uistack(ax1, 'top');

    % Histogram of absolute finite-difference slopes
    % axSlope = axes('Position', [0.82 0.08 0.15 0.30]);
    % hold(axSlope, 'on');
    % 
    % dx = theta1_domain(rightIdx) - theta1_domain(leftIdx);
    % 
    % slopes = zeros(49,1);
    % 
    % for k = 1:999
    %     slopes(k) = ...
    %         (all_likelihoods(k,rightIdx) - ...
    %          all_likelihoods(k,leftIdx)) / dx;
    % end
    % 
    % absSlopes = slopes.^2; % square it instead of abs
    % 
    % % Histogram
    % histogram(axSlope, absSlopes, ...
    %     'NumBins', 15, ...
    %     'FaceColor', curveColor, ...
    %     'FaceAlpha', 0.75, ...
    %     'EdgeColor', 'none');
    % 
    % xlabel(axSlope, '$(\frac{d\log L}{d\theta})^2$', 'Interpreter', 'latex');
    % % ylabel(axSlope, 'Count');
    % 
    % set(axSlope, ...
    %     'FontSize', 10, ...
    %     'LineWidth', 1.2, ...
    %     'TickDir', 'out', ...
    %     'Box', 'off');
    % 
    % xlim(axSlope, [0, 0.03])
    % axSlope.YTick = [];
    % grid(axSlope, 'off');


end


%% MLE FIM relationship - P(thetaMLE given theta*) - Gaussian
figure(12)
clf
hold on;

for j = 1:length(N)
    theta_true = all_thetas(j,:);
    n = N(j);

    r = theta_true(2)/sqrt(n);

    pmle = normpdf(theta1_domain, theta_true(1), r);
    plot(theta1_domain, pmle, 'LineWidth', 3)
end

xline(theta_true(1), 'k--', 'LineWidth', 3)


%% MLE FIM relationship - S(thetaMLE given theta*) - Gaussian
figure(13)
clf
hold on

nSamples = 10000;

for j = 1:length(N)

    theta_true = all_thetas(j,:);
    n = N(j);

    % Theoretical distribution of MLE
    r = theta_true(2)/sqrt(n);

    % Sample MLE directly from theoretical MLE distribution
    mle_samples = theta_true(1) + r*randn(nSamples,1);

    % Score evaluated at the TRUE theta_1
    slopes = n*(mle_samples - theta_true(1)) / theta_true(2)^2;

    % Squared slope
    slopeSquared = slopes.^2;

    % Histogram
    histogram(slopeSquared, ...
        'NumBins', 50, ...
        'Normalization', 'pdf', ...
        'FaceAlpha', 0.7);

    % Empirical mean
    empiricalMean = mean(slopeSquared);

    % Theoretical mean = Fisher information
    theoreticalMean = n / theta_true(2)^2;

    % Plot empirical mean
    xline(empiricalMean, ...
        '--r', ...
        'LineWidth', 2);

    % Plot theoretical mean
    xline(theoreticalMean, ...
        '--k', ...
        'LineWidth', 2);

    % Print values
    fprintf('n = %d, theta_2 = %.2f\n', n, theta_true(2));
    fprintf('Empirical mean  = %.6f\n', empiricalMean);
    fprintf('Theoretical mean = %.6f\n\n', theoreticalMean);

end

xlabel('$(d\log L/d\theta_1)^2$', 'Interpreter', 'latex')
ylabel('Density')

legend('Histogram', ...
       'Empirical mean', ...
       'Theoretical mean', ...
       'Location', 'best')

xlim([0, 0.01])



%% Central Limit Theorem demonstration

% Mixture weights
w = [0.75, 0.25];

% Component means
muComp = [2, 7];

% Component standard deviations
sigmaComp = [0.5, 2];

% Population mean
mu = sum(w .* muComp);

% Population variance
sigma2 = sum(w .* (sigmaComp.^2 + muComp.^2)) - mu^2;
sigma = sqrt(sigma2);


% Simulation settings

nRepeats = 10000;

% Number of observations averaged in each experiment
nValues = [1 5 100];


% Plot the original distribution

x = linspace(-2, 14, 1000);

pdf = w(1) * normpdf(x, muComp(1), sigmaComp(1)) + ...
      w(2) * normpdf(x, muComp(2), sigmaComp(2));

figure(14)
clf
hold on

plot(x, pdf, ...
    'k-', ...
    'LineWidth', 2.5);

xline(mu, ...
    '--r', ...
    'LineWidth', 1.5);

xlabel('$x$', 'Interpreter', 'latex')
ylabel('Density')

title('Starting distribution')

legend('PDF', 'Population mean', ...
    'Location', 'northeast')

box off


% Central Limit Theorem


for j = 1:length(nValues)

    figure(14+j)
    clf

    n = nValues(j);

    % Generate samples

    % Select mixture component for every observation
    component = rand(nRepeats, n) > w(1);

    % Generate observations from the two Gaussian components
    samples = zeros(nRepeats, n);

    idx1 = ~component;
    idx2 = component;

    samples(idx1) = muComp(1) + sigmaComp(1) * randn(sum(idx1(:)), 1);
    samples(idx2) = muComp(2) + sigmaComp(2) * randn(sum(idx2(:)), 1);

    % Calculate mean of each experiment

    sampleMeans = mean(samples, 2);


    % Plot histogram

    nexttile
    hold on

    histogram(sampleMeans, ...
        'NumBins', 50, ...
        'Normalization', 'pdf', ...
        'FaceColor', [0.25 0.45 0.70], ...
        'FaceAlpha', 0.70, ...
        'EdgeColor', 'none');


    % CLT prediction

    sigmaMean = sigma / sqrt(n);

    xTheory = linspace( ...
        min(sampleMeans), ...
        max(sampleMeans), ...
        1000);

    pdfTheory = normpdf( ...
        xTheory, ...
        mu, ...
        sigmaMean);

    plot(xTheory, pdfTheory, ...
        'r-', ...
        'LineWidth', 2);


    % Population mean

    xline(mu, ...
        '--k', ...
        'LineWidth', 1.5);


    % Formatting

    title(sprintf('$n = %d$', n), ...
        'Interpreter', 'latex')

    xlabel('$\bar{X}$', ...
        'Interpreter', 'latex')

    ylabel('Density')

    box off
end



%% Export Figures for Paper Supplimental Figure
outputFolder = 'AnnualReview_Figures';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Overall paper canvas
fullWidth = 6.33;
fullHeight = 7.9;

% 4 x 3 grid
plotWidth = fullWidth / 3;
plotHeight = fullHeight / 3;

for figNum = 10:20

    fig = figure(figNum);

    % Remove figure-level title
    sgtitle(fig, '');

    % Find all axes
    axesList = findall(fig, 'Type', 'Axes');

    for i = 1:length(axesList)

        ax = axesList(i);

        % Remove title
        ax.Title.String = '';
        ax.Title.Visible = 'off';

        % Remove axis labels
        ax.XLabel.String = '';
        ax.XLabel.Visible = 'off';

        ax.YLabel.String = '';
        ax.YLabel.Visible = 'off';

        % Remove ticks and tick labels
        ax.XTick = [];
        ax.YTick = [];

        ax.XTickLabel = [];
        ax.YTickLabel = [];

        % Remove tick marks
        ax.TickLength = [0 0];

    end

    % Remove legends
    legends = findall(fig, 'Type', 'Legend');

    if ~isempty(legends)
        delete(legends);
    end

    % Remove colorbar labels/ticks
    colorbars = findall(fig, 'Type', 'ColorBar');

    for i = 1:length(colorbars)

        cb = colorbars(i);

        cb.TickLabels = [];
        cb.Label.String = '';

    end

    % Set physical dimensions for 4 x 3 grid
    fig.Units = 'inches';
    fig.Position(3:4) = [plotWidth plotHeight];

    % Export
    fileName = sprintf('figure%d.svg', figNum);

    exportgraphics(fig, ...
        fullfile(outputFolder, fileName), ...
        'ContentType', 'vector');

end

disp('Figures 10-20 exported successfully.');




%% MLE FIM relationship - Single cell - Bursting Model
Model_chg.tSpan = linspace(0,25,31);
nCellsInExperiment = zeros(size(Model_chg.tSpan));
nCellsInExperiment([1]) = 1;
Model_chg.parameters{2,2} = star_koff;
Model_chg.parameters{1,2} = final_kon;
Model_chg.parameters{6,2} = star_koff;
Model_chg.parameters{5,2} = final_kon;
Model_chg = Model_chg.solve;
Model_chg.ssaOptions.Nexp = 5000;

Model_chg.fittingOptions.modelVarsToFit = [1];

% Model_chg.plotFSP(plotType='meansAndDevs', SpeciesIdx=[2], Title='testing steady state') % Test successful 
Model_chg.sampleDataFromFSP(saveFile='dataForFIMIntro.csv',nCells=nCellsInExperiment,species2save={'mRNA'});
Model_chg = Model_chg.loadData('dataForFIMIntro.csv', {'mRNA', 'exp1_mRNA'});

kon_domain = logspace(-2,2.5, 200); 
count_domain = 0:200;

pars = [Model_chg.parameters{:,2}];

likelihoods = zeros(size(kon_domain));
for i = 1:length(kon_domain)
    pars(1) = kon_domain(i);
    l = Model_chg.computeLikelihood(pars);
    likelihoods(i) = l;
end

%% MLE FIM relationship - Single cell - Plotting
figure(101);
clf
plot(kon_domain, likelihoods, 'LineWidth', 1.5)
hold on

% ax.Box = 'on';
% ax.LineWidth = 1;

ax = gca;
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.015];

% Find MLE
[~, max_idx] = max(likelihoods);
mle = kon_domain(max_idx);

% True parameter
xline(final_kon, 'k--', 'LineWidth', 2)

% MLE
xline(mle, 'r-', 'LineWidth', 2)

set(gca, 'XScale', 'log')

xlabel('k_{on}')
ylabel('Log-Likelihood')
legend('Likelihood', 'True k_{on}', 'MLE', 'Location', 'best')
% grid on

%% MLE FIM relationship - Multiple cell - Bursting Model - Compute
T = array2table(count_domain, ...
    'VariableNames', compose("exp%d_mRNA", count_domain));

T.time = 0;
T = movevars(T, 'time', 'Before', 1);

writetable(T, 'fakeData.csv');

if false
    likelihoods = zeros([length(count_domain), length(kon_domain)]);
    for i = 1:length(kon_domain)
        pars(1) = kon_domain(i);
        Model_chg.parameters(:,2) = num2cell(pars);
        Model_chg = Model_chg.solve;
        for j = 1:length(count_domain)
            Model_chg = Model_chg.loadData('fakeData.csv', {'mRNA', sprintf('exp%d_mRNA',j-1)});
            l = Model_chg.computeLikelihood(pars, Model_chg.Solutions.stateSpace, false, true);
            likelihoods(j, i) = l;
        end
    end
    save('likeloodFunctions.mat', 'likelihoods')
end


%% MLE FIM relationship - Multiple cell - Bursting Model - Plotting
load("likeloodFunctions.mat")

T = readtable('dataForFIMIntro.csv');
samples = T{1, 2:end};

L = likelihoods(samples(1:20)+1,:);

[~, max_idx] = max(L, [], 2);
mle = kon_domain(max_idx);

figure(102);
clf
hold on

plot(kon_domain, L, 'lineWidth', 1.5)
set(gca, 'XScale', 'log')

ylims = ylim;
ymin = ylims(1);

plot(mle, ymin * ones(size(mle)), 'rx', ...
    'MarkerSize', 12, 'LineWidth', 2)

xline(final_kon, 'k--', 'lineWidth', 2)

ax = gca;
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.015];

%% MLE FIM relationship - MLE Spread
L = likelihoods(samples+1,:);
[~, max_idx] = max(L, [], 2);
mle = kon_domain(max_idx);
mle_log = log(mle);
mleVar_log = var(mle_log);

figure(103);
clf

% Define bins uniformly in log10 space
nBins = 25;
xlims = [min(mle), max(mle)];

log_edges = linspace(log10(xlims(1)), ...
                     log10(xlims(2)), nBins+1);

bin_edges = 10.^log_edges;

histogram(mle, bin_edges, ...
    'Normalization', 'pdf', ...
    'FaceColor', [0.2 0.5 0.8], ...
    'EdgeColor', 'none')

hold on

xline(final_kon, 'k--', 'LineWidth', 2)
xline(mean(mle), 'r', 'LineWidth', 2)

set(gca, 'XScale', 'log')
xlim(xlims)

xlabel('MLE k_{on}')
ylabel('Probability Density')
% grid on

ax = gca;
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.015];

%% MLE FIM relationship - sensitivity and FIM prediction
L = likelihoods(samples+1,:);

% Numerical derivative of each likelihood curve
dL_dkon = zeros(size(L));

for j = 1:size(L,1)
    dL_dkon(j,:) = gradient(L(j,:), kon_domain);
end

% Evaluate sensitivity at the true parameter
sensitivity = zeros(size(samples));

for j = 1:size(L,1)
    sensitivity(j) = interp1(kon_domain, dL_dkon(j,:), ...
                             final_kon, 'linear');
end

% Sensitivity squared
sensitivity_squared = sensitivity.^2;

% Plot histogram
figure(104); 
hold on;
histogram(sensitivity_squared, 50, 'Normalization', 'pdf')
xlabel('Sensitivity^2')
ylabel('Probability Density')
title('Sensitivity Squared at True Parameter')

xlim([0,0.02])

% Empirical MLE variance -> information in log space
% xline(1/mleVar_log, 'b--', 'LineWidth', 2)

% Model FIM
FIM = Model_chg.computeFIM();
FIMEstimate = FIM{1};
xline(FIMEstimate, 'r', 'LineWidth', 2)

xline(mean(sensitivity_squared), 'k--', 'LineWidth', 2)

ax = gca;
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.015];


% legend('Sensitivity^2', ...
%        'Mean sensitivity^2', ...
%        '1 / MLE variance', ...
%        'FIM')



%% MLE and sensitivity for multiple cell numbers

rng(1);  % Reproducibility

% Can contain as many cell numbers as you want
cell_numbers = [2 4 10 50 100 200 500];
nSets = 10000;

% All available single-cell likelihood curves
L_all = likelihoods(samples+1,:);

MLEs = cell(length(cell_numbers),1);
SensitivitySquared = cell(length(cell_numbers),1);

% Calculate MLEs and sensitivity for all cell numbers

for n = 1:length(cell_numbers)

    nCells = cell_numbers(n);

    MLEs{n} = zeros(nSets,1);
    SensitivitySquared{n} = zeros(nSets,1);

    for k = 1:nSets

        % Randomly select cells
        idx = randperm(size(L_all,1), nCells);

        % Sum log-likelihoods across cells
        L_sum = sum(L_all(idx,:), 1);

        % Find MLE
        [~, max_idx] = max(L_sum);
        MLEs{n}(k) = kon_domain(max_idx);

        % Sensitivity of summed log-likelihood
        dL_dkon = gradient(L_sum, kon_domain);

        % Evaluate sensitivity at true kon
        sensitivity = interp1(kon_domain, dL_dkon, ...
                              final_kon, 'linear');

        % Sensitivity squared
        SensitivitySquared{n}(k) = sensitivity^2;
    end

    fprintf('Finished %d cells\n', nCells);
end


% Plot ONLY the first 4 cell numbers

figure(105);
clf

nPlot = min(4, length(cell_numbers));

reference_MLE = MLEs{1};

xlims = [min(reference_MLE), max(reference_MLE)];

nBins = 32;

log_edges = linspace(log10(xlims(1)), ...
                     log10(xlims(2)), nBins+1);

bin_edges = 10.^log_edges;

for n = 1:nPlot

    subplot(2,2,n)

    histogram(MLEs{n}, bin_edges, ...
        'Normalization', 'pdf', ...
        'FaceColor', [0.2 0.5 0.8], ...
        'FaceAlpha', 0.45, ...
        'EdgeColor', 'none')

    hold on

    % True parameter
    xline(final_kon, 'k--', 'LineWidth', 2)

    % Mean MLE
    xline(mean(MLEs{n}), 'r-', 'LineWidth', 2)

    set(gca, 'XScale', 'log')
    xlim(xlims)

    xlabel('MLE k_{on}')
    ylabel('Probability Density')
    title(sprintf('%d Cells', cell_numbers(n)))

    % grid on
end

ax = gca;
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.015];



%% MLE variance and FIM Convergence 
mleVar = zeros(size(cell_numbers));
fimVar = zeros(size(cell_numbers));

for n = 1:length(cell_numbers)

    % Empirical variance of MLE
    mleVar(n) = var(MLEs{n});

    % FIM prediction
    fimVar(n) = 1 / (cell_numbers(n) * FIMEstimate);

end

figure(106);
clf

plot(cell_numbers, mleVar, 'ko-', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'k')
hold on

plot(cell_numbers, fimVar, 'r^-', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'r')

xlabel('Number of Cells')
ylabel('Variance of k_{on}')
legend('MLE variance', 'FIM prediction', ...
    'Location', 'best')

% grid on

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


%% MSE of MLE vs number of cells

mleMSE = zeros(size(cell_numbers));

for n = 1:length(cell_numbers)

    % MLE estimates for this number of cells
    estimates = MLEs{n};

    % Mean squared error relative to true parameter
    mleMSE(n) = mean((estimates - final_kon).^2);

end

% Plot

figure(107);
clf

plot(cell_numbers, mleMSE, 'ko-', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'k')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlabel('Number of Cells')
ylabel('MLE MSE')
title('MLE Mean Squared Error vs Number of Cells')

% grid on

ax = gca;
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.015];



%% MSE between MLE variance and FIM variance

varMSE = (mleVar - fimVar).^2;

figure(108);
clf

plot(cell_numbers, varMSE, 'ko-', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'k')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlabel('Number of Cells')
ylabel('MSE: MLE Variance vs FIM Variance')
title('MSE Between Empirical Variance and FIM Estimate')

% grid on

ax = gca;
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.TickLength = [0.015 0.015];


%% Export Figures for Paper Supplimental Figure
outputFolder = 'AnnualReview_Figures';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Overall paper canvas
fullWidth = 6.33;
fullHeight = 7.9;

% 4 x 3 grid
plotWidth = fullWidth / 3;
plotHeight = fullHeight / 3;

for figNum = 101:108

    fig = figure(figNum);

    % Remove figure-level title
    sgtitle(fig, '');

    % Find all axes
    axesList = findall(fig, 'Type', 'Axes');

    for i = 1:length(axesList)

        ax = axesList(i);

        % Remove title
        ax.Title.String = '';
        ax.Title.Visible = 'off';

        % Remove axis labels
        ax.XLabel.String = '';
        ax.XLabel.Visible = 'off';

        ax.YLabel.String = '';
        ax.YLabel.Visible = 'off';

        % Remove ticks and tick labels
        % ax.XTick = [];
        % ax.YTick = [];

        ax.XTickLabel = [];
        ax.YTickLabel = [];

        % Remove tick marks
        % ax.TickLength = [0 0];

    end

    % Remove legends
    legends = findall(fig, 'Type', 'Legend');

    if ~isempty(legends)
        delete(legends);
    end

    % Remove colorbar labels/ticks
    colorbars = findall(fig, 'Type', 'ColorBar');

    for i = 1:length(colorbars)

        cb = colorbars(i);

        cb.TickLabels = [];
        cb.Label.String = '';

    end

    % Set physical dimensions for 4 x 3 grid
    fig.Units = 'inches';
    fig.Position(3:4) = [plotWidth plotHeight];

    % Export
    fileName = sprintf('figure%d.svg', figNum);

    exportgraphics(fig, ...
        fullfile(outputFolder, fileName), ...
        'ContentType', 'vector');

end

disp('Figures 101-108 exported successfully.');




%%



























% return

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

% Model = Model.solve;

% Model.plotFSP



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


%% PDO - Effect on Distributions
% Pick a parameter set that has an interesting looking PDF.
%                                     PRIOR
Model.parameters = {'kon0',0.01;...  % logn(-1,2)
    'koff0',0.01;...                   % logn(0,2)
    'kr',10;...                     % logn(1,2)
    'g',0.01;...                    % logn(-2,2)
    'kD',10;...                     % logn(1,2)
    'S0',1;...                      % NA (initial input concentration)
    'S1',5};                        % NA (final input concentration)

f1 = figure(1); clf;
Model.fspOptions.bounds = [];
Model = Model.solve(solver='fsp');
Model.plotFSP(figureNums=f1,plotType='marginals',indTimes=length(Model.tSpan),speciesNames='mRNA',Colors={'r'})

% Add a Binomial PDO 
dropOut = 0.9; % fraction dropout
Model_BinomialPDO = Model;
Model_BinomialPDO.pdoOptions.type = 'Binomial';
Model_BinomialPDO.pdoOptions.unobservedSpecies = 'gON';
Model_BinomialPDO.pdoOptions.props.CaptureProbabilityS1 = 0;    % Gene State is not measured
Model_BinomialPDO.pdoOptions.props.CaptureProbabilityS2 = 1-dropOut; % 95% dropout from RNA
[~,Model_BinomialPDO] = Model_BinomialPDO.generatePDO;
hold on
Model_BinomialPDO.plotFSP(figureNums=f1,plotType='marginals',indTimes=length(Model_BinomialPDO.tSpan),...
    speciesNames='mRNA',includePDO=true,Colors={'k'})
% set(gca,'yscale','log','ylim',[1e-5,1])


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

% TODO - make plots of these measurements along the length of input
% TODO - make plots of each optimality vs NCells for each stratagy
% TODO - make plot of FIM-1 for the original experiment


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


%% Plots of FIM predicted uncertainties
freePars = [1:4];
f1 = figure(201);
f2 = figure(202);
f3 = figure(203);

% The following plots the heatmap showing the
Model.plotFIMResults(FIM_Opt^(-1)/log(10)^2, 'log',...
    Model.parameters(freePars,1),...
    [Model.parameters{freePars,2}],...
    PlotEllipses=true,EllipseFigure=f1,...
    FigureHandle=f3,...
    Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]),...
    LogThreshold=-4,...
    HeatMapType='invfim',...
    MatrixType='invfim');


%%
f4 = figure(204)
Model.plotFIMResults(FIM_Opt^(-1)/log(10)^2, 'log',...
    Model.parameters(freePars,1),...
    [Model.parameters{freePars,2}],...
    PlotEllipses=true,EllipseFigure=f4,...
    EllipsePairs=[1,2], ...
    FigureHandle=f3,...
    Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]),...
    LogThreshold=-4,...
    HeatMapType='invfim',...
    MatrixType='invfim');

hold on

C = FIM_Opt^(-1)/log(10)^2;

% Parameters corresponding to your ellipse pair
C2 = C([1 2],[1 2]);

% Eigenvectors/eigenvalues
[V,D] = eig(C2);

% Sort eigenvalues from smallest to largest
[lambda,idx] = sort(diag(D));
V = V(:,idx);

% Center of ellipse
x0 = log10(Model.parameters{2,2});
y0 = log10(Model.parameters{1,2});

% Scale factor for visualization
scale = 2;

% Small eigenvalue direction
quiver(x0,y0,...
    V(2,1)*sqrt(lambda(1))*scale,...
    V(1,1)*sqrt(lambda(1))*scale,...
    0,...
    'LineWidth',2,...
    'Color','r',...
    'MaxHeadSize',0.5);

% Large eigenvalue direction
quiver(x0,y0,...
    V(2,2)*sqrt(lambda(2))*scale,...
    V(1,2)*sqrt(lambda(2))*scale,...
    0,...
    'LineWidth',2,...
    'Color','b',...
    'MaxHeadSize',0.5);

% TODO - Add MLE estimates to plot
% TODO - change exp for this to acheive MLE spread
% TODO - plot FIM-1 for optimatlity descriptions 
%% 
nCellsInExperiment = zeros(size(Model.tSpan));
nCellsInExperiment([1]) = 1;
Model = Model.solve;
Model.ssaOptions.Nexp = 5000;

Model.fittingOptions.modelVarsToFit = [1];

% Model_chg.plotFSP(plotType='meansAndDevs', SpeciesIdx=[2], Title='testing steady state') % Test successful 
Model.sampleDataFromFSP(saveFile='dataForFIMIntro.csv',nCells=nCellsInExperiment,species2save={'mRNA'});
Model = Model.loadData('dataForFIMIntro.csv', {'mRNA', 'exp1_mRNA'});




return
















%%


%%


%%
















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









%% Functions 
function plotHeatmap(M, rowNames, colNames, titleText)

    % =============================================================
    % Signed-logarithmic heatmap
    %
    % ORIGINAL MATRIX M IS NEVER MODIFIED.
    %
    % Negative -> blue
    % Zero     -> white
    % Positive -> red
    %
    % Powers of 10 are equally spaced in colour space.
    % =============================================================

    if size(M,1) ~= numel(rowNames)
        error('Number of row names must equal number of rows.');
    end

    if size(M,2) ~= numel(colNames)
        error('Number of column names must equal number of columns.');
    end

    rowNames = cellstr(rowNames);
    colNames = cellstr(colNames);

    % -------------------------------------------------------------
    % Colours
    % -------------------------------------------------------------

    blue  = [0.05 0.25 0.75];
    white = [1.00 1.00 1.00];
    red   = [0.80 0.10 0.10];

    % -------------------------------------------------------------
    % Get largest magnitude
    % -------------------------------------------------------------

    maxValue = max(abs(M(:)));

    if maxValue == 0
        maxValue = 1;
    end

    % -------------------------------------------------------------
    % Determine decade range
    %
    % Example:
    %
    % maxValue = 4e15
    %
    % gives approximately:
    %
    % 1e12  1e13  1e14  1e15  1e16
    %
    % -------------------------------------------------------------

    maxExponent = ceil(log10(maxValue));

    nDecades = 4;

    minExponent = maxExponent - nDecades;

    % -------------------------------------------------------------
    % Convert M -> COLOR COORDINATE
    %
    % M itself is NOT changed.
    %
    % Coordinate:
    %
    % negative values : [-1,0]
    % zero            : 0
    % positive values : [0,1]
    %
    % The magnitude is logarithmically positioned.
    % -------------------------------------------------------------

    C = zeros(size(M));

    idx = M ~= 0;

    if any(idx(:))

        magnitude = abs(M(idx));

        t = (log10(magnitude) - minExponent) / ...
            (maxExponent - minExponent);

        % Clamp
        t = max(0, min(1, t));

        C(idx) = sign(M(idx)) .* t;

    end

    % -------------------------------------------------------------
    % Convert colour coordinate to RGB
    % -------------------------------------------------------------

    RGB = zeros([size(M), 3]);

    for i = 1:size(M,1)

        for j = 1:size(M,2)

            c = C(i,j);

            if c < 0

                % Blue -> white
                q = abs(c);

                RGB(i,j,:) = ...
                    blue + q .* (white - blue);

            elseif c > 0

                % White -> red
                q = c;

                RGB(i,j,:) = ...
                    white + q .* (red - white);

            else

                % EXACTLY ZERO
                RGB(i,j,:) = white;

            end

        end

    end

    % -------------------------------------------------------------
    % Plot RGB image
    % -------------------------------------------------------------

    image(RGB);

    ax = gca;

    axis image;

    % -------------------------------------------------------------
    % Create a custom colourbar
    %
    % We make a separate invisible image whose colour coordinate
    % runs from -1 to +1.
    % -------------------------------------------------------------

    hold on;

    % Dummy invisible image used only for the colourbar
    dummy = imagesc([-1 1; -1 1]);

    dummy.Visible = 'off';

    % Use the same blue-white-red colormap
    n = 256;

    nBlue = 128;
    nRed  = 128;

    blueMap = [
        linspace(blue(1), white(1), nBlue)', ...
        linspace(blue(2), white(2), nBlue)', ...
        linspace(blue(3), white(3), nBlue)'
    ];

    redMap = [
        linspace(white(1), red(1), nRed)', ...
        linspace(white(2), red(2), nRed)', ...
        linspace(white(3), red(3), nRed)'
    ];

    colormap(ax, [blueMap; redMap]);

    % -------------------------------------------------------------
    % Colorbar
    % -------------------------------------------------------------

    cb = colorbar;

    clim([-1 1]);

    % -------------------------------------------------------------
    % Construct tick positions DIRECTLY.
    %
    % This is the important part:
    %
    % -1, -0.75, -0.5, -0.25, 0, ...
    %
    % are strictly increasing.
    % -------------------------------------------------------------

    exponents = minExponent:maxExponent;

    % Positions corresponding to powers of ten
    %
    % minExponent -> 0
    % maxExponent -> 1

    decadePosition = ...
        (exponents - minExponent) ./ ...
        (maxExponent - minExponent);

    % Negative side
    negativePositions = -fliplr(decadePosition);

    % Positive side
    positivePositions = decadePosition;

    % Combine in increasing order
    tickPositions = [
        negativePositions ...
        0 ...
        positivePositions
    ];

    % Remove duplicate zero if it occurs
    tickPositions = unique(tickPositions, 'sorted');

    % -------------------------------------------------------------
    % Labels
    % -------------------------------------------------------------

    tickLabels = cell(size(tickPositions));

    for k = 1:numel(tickPositions)

        p = tickPositions(k);

        if p == 0

            tickLabels{k} = '0';

        else

            % Recover exponent from position
            e = minExponent + ...
                abs(p) * (maxExponent - minExponent);

            e = round(e);

            if p < 0
                tickLabels{k} = sprintf('$-10^{%d}$', e);
            else
                tickLabels{k} = sprintf('$10^{%d}$', e);
            end

        end

    end

    cb.Ticks = tickPositions;
    cb.TickLabels = tickLabels;

    cb.TickLabelInterpreter = 'latex';
    cb.Label.String = 'Value';

    % -------------------------------------------------------------
    % Axes
    % -------------------------------------------------------------

    ax.XTick = 1:numel(colNames);
    ax.YTick = 1:numel(rowNames);

    ax.XTickLabel = colNames;
    ax.YTickLabel = rowNames;

    ax.FontSize = 11;
    ax.LineWidth = 0.5;
    ax.TickDir = 'out';
    ax.Box = 'on';
    
    hold on
    
    % Vertical cell boundaries
    for x = 0.5:1:size(M,2)+0.5
        plot([x x], [0.5 size(M,1)+0.5], ...
            'Color', [0.5 0.5 0.5], ...
            'LineWidth', 0.5);
    end
    
    % Horizontal cell boundaries
    for y = 0.5:1:size(M,1)+0.5
        plot([0.5 size(M,2)+0.5], [y y], ...
            'Color', [0.5 0.5 0.5], ...
            'LineWidth', 0.5);
    end
    
    hold off

    xlabel('Parameter');
    ylabel('Parameter');

    axis square;

    % -------------------------------------------------------------
    % +/- symbols
    % -------------------------------------------------------------

    for i = 1:size(M,1)

        for j = 1:size(M,2)

            if M(i,j) > 0
                symbol = '+';

            elseif M(i,j) < 0
                symbol = '−';

            else
                symbol = '0';
            end

            text(j, i, symbol, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'Color', 'black');

        end

    end

    % -------------------------------------------------------------
    % Title
    % -------------------------------------------------------------

    if nargin >= 4 && ~isempty(titleText)

        title(titleText, ...
            'FontSize', 14, ...
            'FontWeight', 'normal');

    end

    hold off;

end
