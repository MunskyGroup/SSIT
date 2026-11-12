%% SSIT/Examples/example_15_PGF_STL1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving the CME through its probability generating function (PGF)
%   * Solve the 4-state STL1 yeast model with SSIT.solve(solver='pgf')
%   * Compare the PGF solution with the Finite State Projection (FSP)
%
% The PGF solver represents the distribution of the "large" species by its
% probability generating function
%       G(z,t) = E[ z^mRNA(t) ],
% which obeys a linear PDE in z when the propensities are polynomial in that
% species (mass action).  The PDE is solved by collocation on a Smolyak
% sparse grid of Chebyshev nodes in z, so the cost depends on the number of
% grid nodes, not on the size of the state space.  Small, bounded species
% (here the gene states g1..g4, which satisfy g1+g2+g3+g4 = 1) are detected
% automatically and enumerated: one PGF block is carried per gene state.
%
% The PGF solver returns moments (means, covariances), the probabilities of
% the enumerated states, and G(z,t) itself -- not the full distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% clear
% close all

% Load the pre-set 4-state STL1 model (same as STL1_4state in
% example_01_CreateSSITModels): gene states g1..g4, mRNA, and the
% time-varying Hog1 input that switches on at t0.
STL1_4state = SSIT('STL1_4state');
STL1_4state.summarizeModel

% Times at which the solutions are computed (minutes):
STL1_4state.tSpan = linspace(0,50,101);

% Start from the steady state at t = 0 (before the Hog1 signal).  Both the
% FSP and the PGF solvers follow this option.
STL1_4state.fspOptions.initApproxSS = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Reference solution with the FSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STL1_fsp = STL1_4state;
STL1_fsp.fspOptions.fspTol = 1e-6;     % tight tolerance: used as the reference
STL1_fsp = STL1_fsp.solve(solver='FSP');
tic
STL1_fsp = STL1_fsp.solve(solver='FSP');
fspTime = toc;
fprintf('FSP solve: %.2f s\n', fspTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve the same model with the PGF solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STL1_pgf = STL1_4state;

% Options (see ssit.pgf.defaultOptions for all of them):
%   mu          Smolyak level: more nodes -> more accurate (1D here: 2^mu+1 nodes)
%   lowerBound  z-domain [lowerBound, 1] for mRNA; 'auto' picks it from the
%               solution's mean (re-solving autoPasses times)
%   enumeratedSpecies  'auto' detects g1..g4 from their conservation law
STL1_pgf.pgfOptions.mu = 5;
STL1_pgf.pgfOptions.lowerBound = 'auto';
STL1_pgf.pgfOptions.relTol = 1e-9;
STL1_pgf.pgfOptions.absTol = 1e-12;
STL1_pgf.pgfOptions.verbose = true;

% First solve: translates the model, finds the gene states, chooses the
% z-domain and builds the grid.  These are kept with the model, so later
% solves (e.g. with new parameter values) skip them, as long as the model
% structure is unchanged.
STL1_pgf = STL1_pgf.solve(solver='pgf');
tic; STL1_pgf = STL1_pgf.solve(solver='pgf'); pgfTime = toc;   % timed (warm) solve
pgfSol = STL1_pgf.Solutions.pgf;
% pgfSol.means, pgfSol.covs          : moments of all species vs. time
% pgfSol.enumeratedProbabilities     : P(gene state) vs. time
% pgfSol.Gtot                        : G(z,t) of mRNA at the grid nodes
% pgfSol.grid, pgfSol.lowerBound     : the sparse grid and z-domain used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Compare the PGF and FSP solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = STL1_4state.tSpan;
iR = find(strcmp(STL1_4state.species, 'mRNA'));
geneNames = pgfSol.enumeratedSpecies;
iG = find(ismember(STL1_4state.species, geneNames));

% Moments of the FSP solution (all species, every time):
[fspMeans, fspCovs] = fspMomentsFromSolution(STL1_fsp.Solutions.fsp);
fspSD = sqrt(squeeze(fspCovs(iR,iR,:)))';
pgfSD = sqrt(squeeze(pgfSol.covs(iR,iR,:)))';

% Gene-state probabilities: P(g_i = 1) = E[g_i] in the FSP; the PGF solver
% returns the probability of each enumerated state (columns of
% pgfSol.enumeratedStates), which here are the unit vectors g1..g4.
[~, stateOfGene] = max(pgfSol.enumeratedStates, [], 1);   % state e has gene stateOfGene(e) on
pgfGene = zeros(numel(iG), numel(T));
pgfGene(stateOfGene,:) = pgfSol.enumeratedProbabilities;

% The PGF of mRNA itself, on the grid nodes and on a fine z-grid:
Gfsp_nodes = ssit.pgf.fspToPgf(STL1_fsp.Solutions.fsp, pgfSol.grid.z, iR);
zFine = linspace(pgfSol.lowerBound, 1, 200)';
Gpgf_fine = ssit.pgf.interpolate(pgfSol.grid, pgfSol.Gtot, zFine);
Gfsp_fine = ssit.pgf.fspToPgf(STL1_fsp.Solutions.fsp, zFine, iR);

fprintf('\nPGF vs FSP (%d grid nodes on z in [%.4f, 1], %d gene states)\n', ...
    pgfSol.grid.N, pgfSol.lowerBound, size(pgfSol.enumeratedStates,2));
fprintf('  solve time         : PGF %.2f s,  FSP %.2f s\n', pgfTime, fspTime);
fprintf('  max |G_pgf - G_fsp| at nodes : %.2e\n', max(abs(pgfSol.Gtot - Gfsp_nodes), [], 'all'));
fprintf('  max |G_pgf - G_fsp| on [lower,1] (interpolated): %.2e\n', max(abs(Gpgf_fine - Gfsp_fine), [], 'all'));
fprintf('  max |mean mRNA diff|  : %.2e   (max mean %.1f)\n', ...
    max(abs(pgfSol.means(iR,:) - fspMeans(iR,:))), max(fspMeans(iR,:)));
fprintf('  max |SD mRNA diff|    : %.2e   (max SD %.1f)\n', max(abs(pgfSD - fspSD)), max(fspSD));
fprintf('  max |P(gene) diff|    : %.2e\n', max(abs(pgfGene - fspMeans(iG,:)), [], 'all'));

%% Plots
figure(1); clf
set(gcf, 'Position', [100 100 1300 420]);

% (a) mRNA mean +/- SD
subplot(1,3,1); hold on
h1 = plot(T, fspMeans(iR,:), 'k-', 'LineWidth', 2);
h2 = plot(T, [fspMeans(iR,:) + fspSD; fspMeans(iR,:) - fspSD], 'k--', 'LineWidth', 1);
idx = 1:4:numel(T);
h3 = plot(T(idx), pgfSol.means(iR,idx), 'ro', 'MarkerFaceColor', 'r');
h4 = plot(T(idx), [pgfSol.means(iR,idx) + pgfSD(idx); pgfSol.means(iR,idx) - pgfSD(idx)], 'r^');
xlabel('Time (min)'); ylabel('mRNA count');
title('mRNA mean \pm SD');
legend([h1, h2(1), h3, h4(1)], {'FSP mean','FSP \pm SD','PGF mean','PGF \pm SD'}, 'Location', 'northeast');

% (b) gene-state probabilities
subplot(1,3,2); hold on
cols = lines(numel(iG));
for i = 1:numel(iG)
    plot(T, fspMeans(iG(i),:), '-', 'Color', cols(i,:), 'LineWidth', 2);
end
for i = 1:numel(iG)
    plot(T(idx), pgfGene(i,idx), 'o', 'Color', cols(i,:), 'MarkerFaceColor', cols(i,:));
end
xlabel('Time (min)'); ylabel('Probability');
title('Gene states (lines FSP, markers PGF)');
legend(STL1_4state.species(iG), 'Location', 'east');

% (c) the PGF of mRNA at a few times
subplot(1,3,3); hold on
showT = [1, 11, 21, 41, 101];
cols = parula(numel(showT)+1);
for k = 1:numel(showT)
    it = showT(k);
    plot(zFine, Gfsp_fine(:,it), '-', 'Color', cols(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('t = %g', T(it)));
    plot(pgfSol.grid.z, pgfSol.Gtot(:,it), 'o', 'Color', cols(k,:), ...
        'MarkerFaceColor', cols(k,:), 'HandleVisibility', 'off');
end
xlabel('z'); ylabel('G(z,t) = E[z^{mRNA}]');
title('mRNA PGF (lines FSP, markers PGF nodes)');
legend('Location', 'northwest');

%% Local functions
function [means, covs] = fspMomentsFromSolution(fsp)
% Means and covariances of all species from SSIT FSP solutions
% (Solutions.fsp: one cell per time with an ssit.FspVector in .p).
Nt = numel(fsp);
d = size(fsp{1}.p.data.subs, 2);
means = zeros(d, Nt);
covs = zeros(d, d, Nt);
for it = 1:Nt
    X = fsp{it}.p.data.subs - 1;     % states (subs are 1-based)
    p = fsp{it}.p.data.vals;
    means(:,it) = X' * p;
    covs(:,:,it) = X' * (p .* X) - means(:,it) * means(:,it)';
end
end
