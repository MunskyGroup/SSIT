%% fimExample
% In this script, we show how to set up and solve the FSP-FIM matrix with
% partial observations and probabilistic distortion.
clear all
ModelChoice = 'CentralDogma';  % Two species problem (mRNa and protein)
F2 = SSIT(ModelChoice);
F2.solutionScheme = 'FSP';    % Set solution scheme to FSP.
[FSPsoln,bounds] = F2.solve;  % Solve the FSP analysis
F2.fspOptions.bounds = bounds;% Save bound for faster analyses 
% 
F2.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F2.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

%% Solve for FIM using FSP Sensitivity Results
F2.pdoOptions.unobservedSpecies = [];
F2.pdoOptions.PDO = [];
[fimResults] = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.

cellCounts = 2*ones(size(F2.tSpan));  % Number of cells in each experiment.
[fimTotal,mleCovEstimate,fimMetrics] = F2.evaluateExperiment(fimResults,cellCounts)

%% Compute FIM for partial observations
F2.pdoOptions.PDO=[];
F2.pdoOptions.unobservedSpecies = 'x1';
[fimResults_partialObs] = F2.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[fimTotal_partialObs,mleCovEstimate_partialObs,fimMetrics_partialObs] = F2.evaluateExperiment(fimResults_partialObs,cellCounts);

%% Compute FIM for FIXED Probabilistic Distortion Operator
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

%% THE FOLLOWING ARE NOT YET FUNCTIONAL 

%% Compute FIM for VARIABLE Probabilistic Distortion Operator
% Define probabilistic distortion for each species S1, S2,...
cellCounts = 2*ones(size(F2.tSpan));  % Number of cells in each experiment.
F4 = F2;
F4.pdoOptions.unobservedSpecies = 'x1';
F4.pdoOptions.type = 'Binomial';
F4.pdoOptions.props.CaptureProbabilityS1 = 0.9;
F4.pdoOptions.props.CaptureProbabilityS2 = 0.9;
F4.pdoOptions.props.ParameterGuess = 1; 

% call method to generate the PDO.
F4.pdoOptions.PDO = F4.generatePDO(F4.pdoOptions,[],FSPsoln.fsp,false);
% solve FIM using the specified PDO
[fimResults_BinomialPDO] = F4.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[~,mleCovEstimate,fimMetrics_BinomialPDO] = F4.evaluateExperiment(fimResults_BinomialPDO,cellCounts)
STDs_fixed = sqrt(diag(mleCovEstimate));

% call method to generate the PDO.
F4.pdoOptions.PDO = F4.generatePDO(F4.pdoOptions,[],FSPsoln.fsp,true);
% solve FIM using the specified PDO
[fimResults_BinomialPDO] = F4.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
[~,mleCovEstimate,fimMetrics_BinomialPDO] = F4.evaluateExperiment(fimResults_BinomialPDO,cellCounts)
STDs_free = sqrt(diag(mleCovEstimate));

bar(([STDs_fixed,STDs_free(1:4)]));
% set(gca,'YTick',[0:3],'YTickLabel',10.^[-2:1])



