%% Spatial Compartement FIM Comparison
clear
close all
addpath(genpath('../../../src'));

%% (1) No Spatial Model
%{
          kon        kr     kd
gene_off <-> gene_on -> RNA -> 0
          koff
%}

F1 = SSIT();
F1.parameters = {'kon',10;'koff',1;'kr',1; 'kd', 0.1};
F1.species = {'gene_on';'gene_off'; 'RNA'};
F1.stoichiometry = [1, -1, 0, 0;...
                    -1, 1, 0, 0;...
                    0, 0, 1, -1];
F1.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'kd*RNA'};
F1.initialCondition = [0;1;0;];
F1.summarizeModel
F1 = F1.formPropensitiesGeneral('StandardModel');

% (2) Solve FSP for model
F1.solutionScheme = 'FSP';    % Set solution scheme to FSP.
[FSPsoln,F1.fspOptions.bounds] = F1.solve;  % Solve the FSP analysis
% Plot the results from the FSP analysis
fig11 = figure(11);clf; set(fig11,'Name','Marginal Distributions, gene_on');
fig12 = figure(12);clf; set(fig12,'Name','Marginal Distributions, gene_off');
fig13 = figure(13);clf; set(fig13,'Name','Marginal Distributions, RNA');
F1.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig11,fig12,fig13],{'r-','linewidth',2})

% (3) Solve FSP Sensitivity
F1.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F1.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem
% Plot the results from the sensitivity analysis
fig14 = figure(14);clf; set(fig14,'Name','Marginal Distributions, gene_on');
fig15 = figure(15);clf; set(fig15,'Name','Marginal Distributions, gene_off');
fig16 = figure(16);clf; set(fig16,'Name','Marginal Distributions, RNA');
F1.makePlot(sensSoln,'marginals',[],false,[fig14, fig15, fig16],{'b','linewidth',2})

% (4) Compute FIM using FSP Sensitivity Results
fimResults = F1.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
cellCounts = 10*ones(size(F1.tSpan));  % Number of cells in each experiment.
[fimTotal,mleCovEstimate,fimMetrics] = F1.evaluateExperiment(fimResults,cellCounts)
fig17 = figure(17);clf; set(fig17,'Name','Fim-Predicted Uncertainty Ellipses');
F1.plotMHResults([],fimTotal,'lin',[],fig17)
legend('FIM - Full Observation')

%% (1) Spatial Model
%{
          kon        kr     D               kd
gene_off <-> gene_on -> RNA -> Distance RNA -> 0
          koff
%}

F0 = SSIT();
F0.parameters = {'kon',10;'koff',1;'kr',1; 'D', 1; 'kd', 0.1};
F0.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'};
F0.stoichiometry = [1, -1, 0, 0, 0;...
                    -1, 1, 0, 0, 0;...
                    0, 0, 1, -1, 0;...
                    0, 0, 0, 1, -1;];
F0.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'D*RNA'; 'kd*D_RNA'};
F0.initialCondition = [0;1;0;0;];
F0.summarizeModel
F0 = F0.formPropensitiesGeneral('SpatailModel');

% (2) Solve FSP for model
F0.solutionScheme = 'FSP';    % Set solution scheme to FSP.
[FSPsoln,F0.fspOptions.bounds] = F0.solve;  % Solve the FSP analysis
% Plot the results from the FSP analysis
fig21 = figure(21);clf; set(fig21,'Name','Marginal Distributions, gene_on');
fig22 = figure(22);clf; set(fig22,'Name','Marginal Distributions, gene_off');
fig23 = figure(23);clf; set(fig23,'Name','Marginal Distributions, RNA');
fig24 = figure(24);clf; set(fig24,'Name','Marginal Distributions, Distant RNA');
F0.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig21,fig22,fig23, fig24],{'r-','linewidth',2})

% (3) Solve FSP Sensitivity
F0.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F0.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem
% Plot the results from the sensitivity analysis
fig25 = figure(25);clf; set(fig25,'Name','Marginal Distributions, gene_on');
fig26 = figure(26);clf; set(fig26,'Name','Marginal Distributions, gene_off');
fig27 = figure(27);clf; set(fig27,'Name','Marginal Distributions, RNA');
fig28 = figure(28);clf; set(fig28,'Name','Marginal Distributions, Distant RNA');
F0.makePlot(sensSoln,'marginals',[],false,[fig25,fig26, fig27, fig28],{'b','linewidth',2})

% (4) Compute FIM using FSP Sensitivity Results
fimResults = F0.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
cellCounts = 10*ones(size(F0.tSpan));  % Number of cells in each experiment.
[fimTotal2,mleCovEstimate,fimMetrics] = F0.evaluateExperiment(fimResults,cellCounts)
fig29 = figure(29);clf; set(fig29,'Name','Fim-Predicted Uncertainty Ellipses');
F0.plotMHResults([],fimTotal,'lin',[],fig29)
legend('FIM - Full Observation')