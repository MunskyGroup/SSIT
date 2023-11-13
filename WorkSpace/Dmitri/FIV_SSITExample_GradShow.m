clear all
close all
addpath(genpath('../../src'));

Model = SSIT;
Model.species = {};
Model.initialCondition = [];
Model.propensityFunctions = {};
Model.stoichiometry = [];
Model.parameters = {};
Model.tSpan = linspace(0,2,21);

% Add species

Model = Model.addSpecies('S',100);
Model = Model.addSpecies('Iregr',0);
Model = Model.addSpecies('Iprog',0);
Model = Model.addSpecies('R',0);

Model.stoichiometry = [];

% Add parameters

Model = Model.addParameter({'progTransmit',0.25}); % given effective contact
Model = Model.addParameter({'regrTransmitMultiplier',0.1}); % given effective contact
Model = Model.addParameter({'contact',0.25});
Model = Model.addParameter({'progDeath',1/18});
Model = Model.addParameter({'regrRecovery',1/18});
Model = Model.addParameter({'respawn',1/6});

P = 0.25; % Proportion randomly assigned to progressive and regressive states

% Add reactions

% [S][Iregr][Iprog][R]

% Regressive infection
stoichRegrInfect = [0,-3,1,2];
Model = Model.addReaction('contact*progTransmit*regrTransmitMultiplier*Iregr',stoichRegrInfect);

% Progressive infection
stoichProgInfect = [0,1,-3,2];
Model = Model.addReaction('contact*progTransmit*Iprog',stoichProgInfect);

% Recovery from regressive infection
stoichRegrRecovery = [0,-1,0,1];
Model = Model.addReaction('regrRecovery*Iregr',stoichRegrRecovery);

% Death from progressive infection
stoichProgDeath = [0,0,-1,0];
Model = Model.addReaction('progDeath*Iprog',stoichProgDeath);

% Respawn
stoichRespawn = [1,0,0,0];
Model = Model.addReaction('respawn',stoichRespawn);

Model.summarizeModel
%% Generate Equations for Propensities
%Model = Model.formPropensitiesGeneral('EnzymeModel');

%% Solve the FSP
Model.fspOptions.verbose = true;
[fspSoln,Model.fspOptions.bounds] = Model.solve;

%% Generate In Silico Data
Model.ssaOptions.nSimsPerExpt = 10;
Model.ssaOptions.Nexp = 1;
Model.tSpan = linspace(0,2,9);
Model.sampleDataFromFSP(fspSoln,'InSilicoData_FIV.csv')

%%
Model = Model.loadData('InSilicoData_FIV.csv',{'S','exp1_s4'});
Model.makeFitPlot

%%
% fitOptions = optimset('Display','iter','MaxIter',1000);
% fitOptions.SIG = [];
% for i=1:3
%     fitPars = Model.maximizeLikelihood([],fitOptions);
%     Model.parameters(:,2) = num2cell(fitPars);
% end

%%
% run Metropolis Hastings
MHFitOptions.thin=1;
MHFitOptions.numberOfSamples=1000;
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

% Compute the FIM for each time point
fimResults = Model.computeFIM(sensSoln.sens);

% Choose experiment design (number of measurements per time for each time).
cellCounts = 10*ones(size(Model.tSpan));

% Compute the total FIM:
[fimTotal,mleCovEstimate,fimMetrics] = ...
    Model.evaluateExperiment(fimResults,cellCounts);

% Compare to the MH results
Model.plotMHResults(MHResults,fimTotal)

% Optimize to find the best experiment for same number of spot/cells
allowableCellNumber = 500;
OptimumExperiment = Model.optimizeCellCounts(fimResults,allowableCellNumber,...
    'Determinant',[],[],[]);

% Compute the total FIM:
fimTotalOptimized = Model.evaluateExperiment(fimResults,OptimumExperiment);

% Compare to the MH results
Model.plotMHResults(MHResults,[fimTotal,fimTotalOptimized])

%% trouble shooting unidentifiable models
[Vecs,Vals] = eig(full(fimTotal{1}));
JUncertain = find(diag(Vals)<1e-9)
figure
bar(Vecs(:,JUncertain))