clear all
addpath(genpath('../../src'));

Model = SSIT;
Model.species = {};
Model.initialCondition = [];
Model.propensityFunctions = {};
Model.propensityFunctions = {};
Model.stoichiometry = [];
Model.parameters = {};
Model.tSpan = linspace(0,2,21);

% Add species
for i=1:5
    if i==1
        Model = Model.addSpecies(['E',num2str(i)],1);
    else
        Model = Model.addSpecies(['E',num2str(i)],0);
    end
end

Model = Model.addSpecies('G',0);

Model.stoichiometry = [];

Model = Model.addParameter({'alph',0.1;...
    'bet',0.2});
  % Add reactions
for i = 1:4
    % Model = Model.addParameter({['alpha',num2str(i)],0.1;...
    %                             ['beta',num2str(i+1)],0.2;...
    %                             ['k',num2str(i)],10});
    Model = Model.addParameter({['k',num2str(i)],10});

    stoichForward = zeros(size(Model.species,1),1);
    stoichForward([i,i+1])=[-1,1];
    % Model = Model.addReaction(['alpha',num2str(i),'*E',num2str(i)],stoichForward);
    Model = Model.addReaction(['alph*E',num2str(i)],stoichForward);
    
    stoichBackward = zeros(size(Model.species,1),1);
    stoichBackward([i,i+1])=[1,-1];
    Model = Model.addReaction(['bet*E',num2str(i+1),'*',num2str(i)],stoichBackward);  

    stoichProduction = zeros(size(Model.species,1),1);
    stoichProduction(end)=1;
    Model = Model.addReaction(['k',num2str(i),'*E',num2str(i)],stoichProduction);

end
Model = Model.addParameter({'k5',0.1;...
    'deg',0.2});

stoichProduction = zeros(size(Model.species,1),1);
stoichProduction(end)=1;
Model = Model.addReaction('k5*E5',stoichProduction);

stoichDegradation = zeros(size(Model.species,1),1);
stoichDegradation(end)=-1;
Model = Model.addReaction('deg*G',stoichDegradation);

Model.summarizeModel
%% Generate Equations for Propensities
Model = Model.formPropensitiesGeneral('EnzymeModel');

%% Solve the FSP
Model.fspOptions.verbose = true;
[fspSoln,Model.fspOptions.bounds] = Model.solve;

%% Generate In Silico Data
Model.ssaOptions.nSimsPerExpt = 10;
Model.ssaOptions.Nexp = 1;
Model.tSpan = linspace(0,2,9);
Model.sampleDataFromFSP(fspSoln,'InSilicoData.csv')

%%
Model = Model.loadData('InSilicoData.csv',{'G','exp1_s6'});
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
Model.fittingOptions.modelVarsToFit = [1:8];
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

% Set model to ignore unobservable species.
Model.pdoOptions.unobservedSpecies = {'E1','E2','E3','E4','E5'};
% Model.pdoOptions.unobservedSpecies = {'E1','E2','E3','E4'};

% Compute the FIM for each time point
fimResults = Model.computeFIM(sensSoln.sens);

% Choose experiment design (number of measurements per time for each time).
cellCounts = 100*ones(size(Model.tSpan));

% Compute the total FIM:
[fimTotal,mleCovEstimate,fimMetrics] = ...
    Model.evaluateExperiment(fimResults,cellCounts);

% Compare to the MH results
Model.plotMHResults(MHResults,fimTotal)

% Optimize to find the best experiment for same number of spot/cells
allowableCellNumber = 2100;
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



