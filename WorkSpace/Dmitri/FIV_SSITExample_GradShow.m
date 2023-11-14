clear all
close all
addpath(genpath('../../src'));

Model = SSIT;

% Model.species = {'S','RegrI','ProgI','R'};
% Model.initialCondition = [180;5;5;10];
Model.species = {};
Model.initialCondition = [];

Model.propensityFunctions = {};
Model.stoichiometry = [];
Model.parameters = {};
Model.tSpan = linspace(0,20,4);

% Add species

Model = Model.addSpecies('S',20);
Model = Model.addSpecies('RegrI',4);
Model = Model.addSpecies('ProgI',4);
Model = Model.addSpecies('R',8);

Model.stoichiometry = [];

% Add parameters

Model = Model.addParameter({'progTransmit',0.25}); % given effective contact
Model = Model.addParameter({'regrTransmitMultiplier',0.1}); % given effective contact
Model = Model.addParameter({'contact',0.25});
Model = Model.addParameter({'progDeath',1/18});
Model = Model.addParameter({'regrRecovery',1/18});
Model = Model.addParameter({'respawn',1/6});
% Model = Model.addParameter({'infectionArrival',1/10});

P = 0.25; % Proportion randomly assigned to progressive and regressive states

% Add reactions

% [S][RegrI][ProgI][R]
% 
% % Regressive infection due to arrival from outside
% stoichRegrInfectOutside = [-1;1;0;0];
% Model = Model.addReaction('infectionArrival*S/4',stoichRegrInfectOutside);
% 
% % Progressive infection due to arrival from outside
% stoichProgInfectOutside = [-1;0;1;0];
% Model = Model.addReaction('infectionArrival*S/4',stoichProgInfectOutside);
% 
% % Abortive infection due to arrival from outside
% stoichAborInfectOutside = [-1;0;0;1];
% Model = Model.addReaction('infectionArrival*S/2',stoichAborInfectOutside);

% Regressive infection from regressive
stoichRegrInfectRegr = [-1;1;0;0];
Model = Model.addReaction('contact*progTransmit*regrTransmitMultiplier*S*RegrI/4',stoichRegrInfectRegr);

% Progressive infection from regressive
stoichProgInfectRegr = [-1;0;1;0];
Model = Model.addReaction('contact*progTransmit*regrTransmitMultiplier*S*RegrI/4',stoichProgInfectRegr);

% Abortive infection from regressive
stoichAborInfectRegr = [-1;0;0;1];
Model = Model.addReaction('contact*progTransmit*regrTransmitMultiplier*S*RegrI/2',stoichAborInfectRegr);

% Regressive infection from progressive
stoichRegrInfectProg = [-1;1;0;0];
Model = Model.addReaction('contact*progTransmit*S*ProgI/4',stoichRegrInfectProg);

% Progressive infection from progressive
stoichProgInfectProg = [-1;0;1;0];
Model = Model.addReaction('contact*progTransmit*S*ProgI/4',stoichProgInfectProg);

% Abortive infection from progressive
stoichAborInfectProg = [-1;0;0;1];
Model = Model.addReaction('contact*progTransmit*S*ProgI/2',stoichAborInfectProg);

% Recovery from regressive infection
stoichRegrRecovery = [0;-1;0;1];
Model = Model.addReaction('regrRecovery*RegrI',stoichRegrRecovery);

% Death from progressive infection
stoichProgDeath = [0;0;-1;0];
Model = Model.addReaction('progDeath*ProgI',stoichProgDeath);

% Respawn
stoichRespawn = [1;0;0;0];
Model = Model.addReaction('respawn',stoichRespawn);

Model.summarizeModel
%% Generate Equations for Propensities
Model = Model.formPropensitiesGeneral('DiseaseSpreadModel');

%% Solve the FSP

Model.tSpan = linspace(0,20,4);
Model.fspOptions.verbose = true;
[fspSoln,Model.fspOptions.bounds] = Model.solve;

%% Generate In Silico Data
Model.ssaOptions.nSimsPerExpt = 50;
Model.ssaOptions.Nexp = 1;
Model.sampleDataFromFSP(fspSoln,'InSilicoData.csv')

%%
Model = Model.loadData('InSilicoData.csv',{'R','exp1_s4'});
Model.makeFitPlot

%%
fitOptions = optimset('Display','iter','MaxIter',1000);
fitOptions.SIG = [];
for i=1:3
    fitPars = Model.maximizeLikelihood([],fitOptions);
    Model.parameters(:,2) = num2cell(fitPars);
end

%%
% run Metropolis Hastings
MHFitOptions.thin=1;
MHFitOptions.numberOfSamples=2000;
MHFitOptions.burnIn=0;
MHFitOptions.progress=true;
MHFitOptions.useFIMforMetHast =true;
MHFitOptions.suppressFSPExpansion = true;
MHFitOptions.CovFIMscale = 1;
MHFitOptions.numChains = 1;
MHFitOptions.saveFile = 'TMPMHChain.mat';
Model.fittingOptions.modelVarsToFit = [1:5];
delete('TMPMHChain.mat')
[newPars,~,MHResults] = Model.maximizeLikelihood(...
    [], MHFitOptions, 'MetropolisHastings');
Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(newPars);
delete('TMPMHChain.mat')

%%
Model.plotMHResults(MHResults)

%%  FIM Calculation
%%       Solve FSP Sensitivity
Model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = Model.solve(fspSoln.stateSpace);  % Solve the sensitivity problem

% All species are observable

% Compute the FIM for each time point
fimResults = Model.computeFIM(sensSoln.sens);

% Choose experiment design (number of measurements per time for each time).
cellCounts = Model.dataSet.nCells;

% Compute the total FIM:
[fimTotal,mleCovEstimate,fimMetrics] = ...
    Model.evaluateExperiment(fimResults,cellCounts);

% Compare to the MH results
Model.plotMHResults(MHResults,fimTotal)

%% Optimize to find the best experiment for same number of spot/cells
allowableCellNumber = sum(Model.dataSet.nCells);
OptimumExperiment = Model.optimizeCellCounts(fimResults,allowableCellNumber,...
    'Determinant',[],[],[])

% Compute the total FIM:
fimTotalOptimized = Model.evaluateExperiment(fimResults,OptimumExperiment);

% Compare to the MH results
Model.plotMHResults(MHResults,[fimTotal,fimTotalOptimized])

%% trouble shooting unidentifiable models
[Vecs,Vals] = eig(full(fimTotal{1}));
JUncertain = find(diag(Vals)<1e-9)
figure
bar(Vecs(:,JUncertain))

%%
figure
subplot(1,2,1)
pie(Model.dataSet.nCells)
subplot(1,2,2)
pie(OptimumExperiment)

%%
det(fimTotalOptimized{1}^-1)
det(fimTotal{1}^-1)