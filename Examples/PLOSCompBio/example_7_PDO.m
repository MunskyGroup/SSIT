%% example_7_PDO
% Example script to handle distorted data with a binomial conditional
% probability distribution operator (PDO).
addpath(genpath('../src'));

%% Preliminaries
% Load our models from example_1_CreateSSITModels; 
% Compute FSP solutions using example_2_SolveSSITModels_FSP and 
% FIMs using example_4_FIM 
%% Comment out the following 3 lines if example_4_FIM has already been run:
clear
close all
example_4_FIM

% View model summaries
Model_FIM.summarizeModel
STL1_FIM.summarizeModel

% Make copies of Model and STL1 Model
Model_PDO = Model_FIM;
STL1_PDO = STL1_FIM;

%% STEP4 == Define Binomial Probabilistic Distortion Operator
Model_PDO.pdoOptions.type = 'Binomial';
Model_PDO.pdoOptions.props.CaptureProbabilityS1 = 0;  % Distortion for 'offGene' (unobserved)
Model_PDO.pdoOptions.props.CaptureProbabilityS2 = 0;  % Distortion for 'offGene' (unobserved)
Model_PDO.pdoOptions.props.CaptureProbabilityS3 = 0.7;% Distortion for 'mRNA' (observed)
Model_PDO.pdoOptions.PDO = Model_PDO.generatePDO(Model_PDO.pdoOptions,[],Model_sensSoln.sens.data,true);
figure(20); contourf(log10(Model_PDO.pdoOptions.PDO.conditionalPmfs{3}),30); colorbar
xlabel('"true" number of mRNA'); ylabel('observed number of mRNA'); set(gca,'fontsize',15);

%% STEP5 == Apply PDO to FSP and Sensitivity Calculations
Model_PDO.solutionScheme = 'FSP'; % Set solution scheme to FSP.
Model_PDO.makePlot(Model_FSPsoln,'marginals',[1:100:301],true,[1,2,3],{'linewidth',2})  % Plot Distorted Marginals
Model_PDO.solutionScheme = 'fspSens'; % Set solution scheme to Sensitivity
Model_PDO.makePlot(Model_sensSoln,'marginals',[1:100:301],true,[3+(1:12)],{'linewidth',2})    % Plot Distorted Sensitivities

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