%% example_5_FIMCalculation
% Example script to set up and solve the FSP-FIM matrix  
% with partial observations and probabilistic distortion.
clear
close all
addpath(genpath('../src'));

%% Preliminaries
% Load our models described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
example_2_SolveSSITModels_FSP
Model.summarizeModel
STL1Model.summarizeModel

% Make copies of Model and STL1 Model
Model_sens = Model;
STL1_sens = STL1Model;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity
Model_sens.solutionScheme = 'fspSens'; 
STL1_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem
[sensSoln,bounds] = Model_sens.solve(FSPsoln.stateSpace); 
[STL1_sensSoln,bounds] = STL1_sens.solve(STL1_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis
% Model:
fig1 = figure(1);clf; set(fig1,'Name','Marginal Sensitivity, offGene');
fig2 = figure(2);clf; set(fig2,'Name','Marginal Sensitivity, onGene');
fig3 = figure(3);clf; set(fig3,'Name','Marginal Sensitivity, mRNA');
Model_sens.makePlot(sensSoln,'marginals',[],false,...
                    [fig1,fig2,fig3],{'b','linewidth',2})
% STL1 Model:
fig4 = figure(1);clf; set(fig4,'Name','Marginal Sensitivity, offGene');
fig5 = figure(2);clf; set(fig5,'Name','Marginal Sensitivity, onGene');
fig6 = figure(3);clf; set(fig6,'Name','Marginal Sensitivity, mRNA');
STL1_sens.makePlot(STL1_sensSoln,'marginals',[],false,...
                   [fig4,fig5,fig6],{'b','linewidth',2})

%% Compute FIMs using FSP sensitivity results
% Model:
% Compute the FIM for full observations
Model_FIM = Model_sens;
fimResults = Model_FIM.computeFIM(sensSoln.sens); 
% cellCounts = 10*ones(size(F2.tSpan));  % Number of cells in each experiment.
% [fimTotal,mleCovEstimate,fimMetrics] = F2.evaluateExperiment(fimResults,cellCounts)
% fig5 = figure(5);clf; set(fig5,'Name','Fim-Predicted Uncertainty Ellipses');
% F2.plotMHResults([],fimTotal,'lin',[],fig5)
% legend('FIM - Full Observation')

% STL1 Model:
% Compute the FIM for full observations
STL1_FIM = STL1_sens;
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); % Compute the FIM for full observations and no distortion.
% cellCounts = 10*ones(size(F2.tSpan));  % Number of cells in each experiment.
% [fimTotal,mleCovEstimate,fimMetrics] = F2.evaluateExperiment(fimResults,cellCounts)
% fig5 = figure(5);clf; set(fig5,'Name','Fim-Predicted Uncertainty Ellipses');
% F2.plotMHResults([],fimTotal,'lin',[],fig5)
% legend('FIM - Full Observation')

%% (5) Compute FIM for Partial Observations
F2.pdoOptions.PDO=[];
F2.pdoOptions.unobservedSpecies = 'x1';
[fimResults_partialObs] = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[fimTotal_partialObs,mleCovEstimate_partialObs,fimMetrics_partialObs] = F2.evaluateExperiment(fimResults_partialObs,cellCounts)
fig6 = figure(6);clf; set(fig6,'Name','Fim-Predicted Uncertainty Ellipses');
F2.plotMHResults([],[fimTotal,fimTotal_partialObs],'lin',[],fig6)
legend('FIM - Full Observation','FIM - Protein Only')

%% (6) Compute FIM for Distorted Observation (Probabilistic Distortion Operator)
F2.pdoOptions.unobservedSpecies = 'x1';
pdoOptions.type = 'Binomial';
% Need to define loss parameter for each species S1, S2,...
pdoOptions.props.CaptureProbabilityS1 = 0;  % Use zero for unobserved species.
pdoOptions.props.CaptureProbabilityS2 = 0.9;

% Call method to generate the PDO.
F2.pdoOptions.PDO = F2.generatePDO(pdoOptions,[],FSPsoln.fsp);

% Plot the PDO
N = size(F2.pdoOptions.PDO.conditionalPmfs{1});
fig7 = figure(7); set(fig7,'Name','Probabilistic Distortion Operator for Protein');
contourf([0:N(1)-1],[0:N(2)-1],log10(F2.pdoOptions.PDO.conditionalPmfs{1}));
% Here we wanted the first PDO for 'x2' because 'x1' was unobserved.
xlabel('Actual');ylabel('Observable');colorbar;

% Solve FIM using the specified PDO
[fimResults_BinomialPDO] = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[fimTotal_BinomialPDO,mleCovEstimate_BinomialPDO,fimMetrics_BinomialPDO] =...
    F2.evaluateExperiment(fimResults_BinomialPDO,cellCounts)
fig8 = figure(8);clf; set(fig8,'Name','Fim-Predicted Uncertainty Ellipses');
F2.plotMHResults([],[fimTotal,fimTotal_partialObs,fimTotal_BinomialPDO],'lin',[],fig8)
legend('FIM - Full Observation','FIM - Protein Only','FIM - Protein with Error')
