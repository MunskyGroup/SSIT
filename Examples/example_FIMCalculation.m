%% fimExample
% In this script, we show how to set up and solve the FSP-FIM matrix with
% partial observations and probabilistic distortion.
clear all

%% (1) Set up Model
ModelChoice = 'CentralDogma';  % Two species problem (mRNa and protein)
F2 = SSIT(ModelChoice);
F2 = F2.formPropensitiesGeneral('FIMExample');

%% (2) Solve FSP for model
F2.solutionScheme = 'FSP';    % Set solution scheme to FSP.
[FSPsoln,F2.fspOptions.bounds] = F2.solve;  % Solve the FSP analysis

%% (3) Solve FSP Sensitivity
F2.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F2.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

%% (4) Compute FIM using FSP Sensitivity Results
fimResults = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
cellCounts = 10*ones(size(F2.tSpan));  % Number of cells in each experiment.
[fimTotal,mleCovEstimate,fimMetrics] = F2.evaluateExperiment(fimResults,cellCounts)

%% (5) Compute FIM for Partial Observations
F2.pdoOptions.PDO=[];
F2.pdoOptions.unobservedSpecies = 'x1';
[fimResults_partialObs] = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[fimTotal_partialObs,mleCovEstimate_partialObs,fimMetrics_partialObs] = F2.evaluateExperiment(fimResults_partialObs,cellCounts)

%% (6) Compute FIM for Distorted Observation (Probabilistic Distortion Operator)
F2.pdoOptions.unobservedSpecies = 'x1';
pdoOptions.type = 'Binomial';
% Need to define loss parameter for each species S1, S2,...
pdoOptions.props.CaptureProbabilityS1 = 0;  % Use zero for unobserved species.
pdoOptions.props.CaptureProbabilityS2 = 0.9;

% call method to generate the PDO.
F2.pdoOptions.PDO = F2.generatePDO(pdoOptions,[],FSPsoln.fsp);

% plot the PDO
N = size(F2.pdoOptions.PDO.conditionalPmfs{1});
figure;
contourf([0:N(1)-1],[0:N(2)-1],log10(F2.pdoOptions.PDO.conditionalPmfs{1}));
% Here we wanted the first PDO for 'x2' because 'x1' was unobserved.
xlabel('Actual');ylabel('Observable');colorbar;

% solve FIM using the specified PDO
[fimResults_BinomialPDO] = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[fimTotal_BinomialPDO,mleCovEstimate_BinomialPDO,fimMetrics_BinomialPDO] =...
    F2.evaluateExperiment(fimResults_BinomialPDO,cellCounts)