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
[STL1_sensSoln,STL1_bounds] = STL1_sens.solve(STL1_FSPsoln.stateSpace); 

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
cellCounts = 10*ones(size(Model_FIM.tSpan));  % Number of cells in each experiment.
[fimTotal,mleCovEstimate,fimMetrics] = Model_FIM.evaluateExperiment(fimResults,cellCounts)
fig7 = figure(7);clf; set(fig7,'Name','Fim-Predicted Uncertainty Ellipses');
Model_FIM.plotMHResults([],fimTotal,'lin',[],fig7)
legend('FIM - Full Observation')

% STL1 Model:
% Compute the FIM for full observations
STL1_FIM = STL1_sens;
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); % Compute the FIM for full observations and no distortion.
STL1_cellCounts = 10*ones(size(STL1_FIM.tSpan));  % Number of cells in each experiment.
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)
fig8 = figure(5);clf; set(fig8,'Name','Fim-Predicted Uncertainty Ellipses');
STL1_FIM.plotMHResults([],STL1_fimTotal,'lin',[],fig8)
legend('FIM - Full Observation')

%% Compute FIM for partial observations
% Model:
Model_PDO = Model_FIM;
Model_PDO.pdoOptions.PDO=[];
Model_PDO.pdoOptions.unobservedSpecies = 'offGene';
[fimResults_partialObs] = Model_PDO.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[fimTotal_partialObs,mleCovEstimate_partialObs,fimMetrics_partialObs] = Model_PDO.evaluateExperiment(fimResults_partialObs,cellCounts)
fig9 = figure(9);clf; set(fig9,'Name','Fim-Predicted Uncertainty Ellipses');
Model_PDO.plotMHResults([],[fimTotal,fimTotal_partialObs],'lin',[],fig6)
legend('FIM - Full Observation','FIM - Protein Only')

% STL1 Model:
STL1_PDO = STL1_FIM;
STL1_PDO.pdoOptions.PDO=[];
STL1_PDO.pdoOptions.unobservedSpecies = 'offGene';
[STL1_PDO_fimResults_partialObs] = STL1_PDO.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[STL1_fimTotal_partialObs,STL1_mleCovEstimate_partialObs,STL1_fimMetrics_partialObs] = STL1_PDO.evaluateExperiment(STL1_PDO_fimResults_partialObs,STL1_cellCounts)
fig10 = figure(10);clf; set(fig10,'Name','Fim-Predicted Uncertainty Ellipses');
STL1_PDO.plotMHResults([],[STL1_fimTotal,STL1_fimTotal_partialObs],'lin',[],fig10)
legend('FIM - Full Observation','FIM - Protein Only')

%% Compute FIM for distorted observation (Probabilistic Distortion Operator)
%% Model:
Model_PDO.pdoOptions.unobservedSpecies = 'offGene';

% Choose the conditional distribution (PDO)
pdoOptions.type = 'Binomial';

% Need to define loss parameter for each species S1, S2,...
pdoOptions.props.CaptureProbabilityS1 = 0;  % Use zero for unobserved species.
pdoOptions.props.CaptureProbabilityS2 = 0.9;

% Call method to generate the PDO
Model_PDO.pdoOptions.PDO = Model_PDO.generatePDO(pdoOptions,[],FSPsoln.fsp);

% Plot the PDO
N = size(Model_PDO.pdoOptions.PDO.conditionalPmfs{1});
fig11 = figure(11); set(fig11,'Name','Probabilistic Distortion Operator for Protein');
contourf([0:N(1)-1],[0:N(2)-1],log10(Model_PDO.pdoOptions.PDO.conditionalPmfs{1}));

% Here we wanted the first PDO for 'x2' because 'x1' was unobserved.
xlabel('Actual');ylabel('Observable');colorbar;

% Solve FIM using the specified PDO
[fimResults_BinomialPDO] = Model_PDO.computeFIM(sensSoln.sens); 
[fimTotal_BinomialPDO,mleCovEstimate_BinomialPDO,fimMetrics_BinomialPDO] =...
    Model_PDO.evaluateExperiment(fimResults_BinomialPDO,cellCounts)
fig12 = figure(12);clf; set(fig12,'Name','Fim-Predicted Uncertainty Ellipses');
Model_PDO.plotMHResults([],[fimTotal,fimTotal_partialObs,fimTotal_BinomialPDO],'lin',[],fig12)
legend('FIM - Full Observation','FIM - Protein Only','FIM - Protein with Error')

%% STL1 Model:
STL1_Model_PDO.pdoOptions.unobservedSpecies = 'offGene';

% Choose the conditional distribution (PDO)
pdoOptions.type = 'Binomial';

% Need to define loss parameter for each species S1, S2,...
pdoOptions.props.CaptureProbabilityS1 = 0;  % Use zero for unobserved species.
pdoOptions.props.CaptureProbabilityS2 = 0.9;

% Call method to generate the PDO
STL1_Model_PDO.pdoOptions.PDO = STL1_Model_PDO.generatePDO(pdoOptions,[],STL1_FSPsoln.fsp);

% Plot the PDO
STL1_N = size(STL1_Model_PDO.pdoOptions.PDO.conditionalPmfs{1});
fig13 = figure(13); set(fig13,'Name','Probabilistic Distortion Operator for Protein');
contourf([0:N(1)-1],[0:N(2)-1],log10(STL1_Model_PDO.pdoOptions.PDO.conditionalPmfs{1}));

% Here we wanted the first PDO for 'x2' because 'x1' was unobserved.
xlabel('Actual');ylabel('Observable');colorbar;

% Solve FIM using the specified PDO
[STL1_fimResults_BinomialPDO] = STL1_Model_PDO.computeFIM(STL1_sensSoln.sens); 
[STL1_fimTotal_BinomialPDO,STL1_mleCovEstimate_BinomialPDO,STL1_fimMetrics_BinomialPDO] =...
    STL1_Model_PDO.evaluateExperiment(STL1_fimResults_BinomialPDO,STL1_cellCounts)
fig14 = figure(14);clf; set(fig14,'Name','Fim-Predicted Uncertainty Ellipses');
STL1_Model_PDO.plotMHResults([],[STL1_fimTotal,STL1_fimTotal_partialObs,STL1_fimTotal_BinomialPDO],'lin',[],fig14)
legend('FIM - Full Observation','FIM - Protein Only','FIM - Protein with Error')