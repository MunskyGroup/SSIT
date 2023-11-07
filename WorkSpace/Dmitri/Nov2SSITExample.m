clear all
addpath(genpath('../../src'));

Model = SSIT;
Model.species = {'E1'};
Model.initialCondition = [1];
Model.propensityFunctions = {};
Model.propensityFunctions = {'kb'}
Model.stoichiometry = [1];
Model.parameters = ({'kb',0.0});
Model.tSpan = linspace(0,2,21);

% Add species
for i=2:5
    Model = Model.addSpecies(['E',num2str(i)],0);
end

Model = Model.addSpecies('G',0);

% Add reactions
for i = 1:4
    Model = Model.addParameter({['alpha',num2str(i)],0.1;...
                                ['beta',num2str(i+1)],0.2;...
                                ['k',num2str(i)],10});

    stoichForward = zeros(size(Model.species,1),1);
    stoichForward([i,i+1])=[-1,1];
    Model = Model.addReaction(['alpha',num2str(i),'*E',num2str(i)],stoichForward);
    
    stoichBackward = zeros(size(Model.species,1),1);
    stoichBackward([i,i+1])=[1,-1];
    Model = Model.addReaction(['beta',num2str(i+1),'*E',num2str(i+1)],stoichBackward);  

    stoichProduction = zeros(size(Model.species,1),1);
    stoichProduction(end)=1;
    Model = Model.addReaction(['k',num2str(i),'*E',num2str(i)],stoichProduction);

end
Model = Model.addParameter({'k5',0.1;...
    'deg',0.2});

stoichProduction = zeros(size(Model.species,1),1);
stoichProduction(end)=1;
Model = Model.addReaction('k5*E5',stoichProduction);

stoichDegradtion = zeros(size(Model.species,1),1);
stoichDegradtion(end)=-1;
Model = Model.addReaction('deg*G',stoichDegradtion);

Model.summarizeModel
%% Generate Equations for Propensities
Model = Model.formPropensitiesGeneral('EnzymeModel');

%% Solve the FSP
Model.fspOptions.verbose = true;
[fspSoln,Model.fspOptions.bounds] = Model.solve;

%% Generate Fake Data
Model.ssaOptions.nSimsPerExpt = 10;
Model.ssaOptions.Nexp = 1;
Model.sampleDataFromFSP(fspSoln,'fakeData.csv')

%%
Model = Model.loadData('fakeData.csv',{'G','exp1_s6'});
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
MHFitOptions.numberOfSamples=1000;
MHFitOptions.burnIn=0;
MHFitOptions.progress=true;
MHFitOptions.useFIMforMetHast =true;
MHFitOptions.CovFIMscale = 1.0;
MHFitOptions.numChains = 1;
MHFitOptions.saveFile = 'TMPMHChain.mat';
Model.fittingOptions.modelVarsToFit = [1,2];
[newPars,~,MHResults] = Model.maximizeLikelihood(...
    [], MHFitOptions, 'MetropolisHastings');
Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(newPars);
delete('TMPMHChain.mat')

%%  FIM Calculation
Model.unobservedSpecies = {'x1',;


