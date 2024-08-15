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

%% Define limits of nucleus size, reaction order of secondary nucleation:
MIN_NUCLEUS_SIZE = 3;
assert(MIN_NUCLEUS_SIZE > 1, ...
    'Minimum nucleus size must be more than a monomer!')
MAX_NUCLEUS_SIZE = 10;
assert(MAX_NUCLEUS_SIZE >= MIN_NUCLEUS_SIZE, ...
    'Maximum nucleus size must be no smaller than minimum size!')
nc = MIN_NUCLEUS_SIZE;
n2 = 2;
assert(n2 > 0, 'Secondary nucleation must proceed via a nucleus!')

%% Define initial monomer population:
INITIAL_MONOMER_POPULATION = 100;
assert(INITIAL_MONOMER_POPULATION >= MIN_NUCLEUS_SIZE, ...
    'There must be sufficiently many monomers to form one nucleus')

%% Add species: xk denotes a k-mer, from 1 up to the maximum considered.

for i = 1 : MAX_NUCLEUS_SIZE
    if i == 1
        Model = Model.addSpecies(['x', num2str(i)], ...
            INITIAL_MONOMER_POPULATION);
    else
        Model = Model.addSpecies(['x', num2str(i)], 0);
    end
end

%% Add parameters

Model.stoichiometry = [];

Model = Model.addParameter({ ...
    'kn', 0.1; ... % Primary nucleation
    'k2', 0.2; ... % Surface-catalyzed secondary nucleation
    'kplus', 0.3; ... % Elongation (at end of linear chain)
    'kminus',0.1 ... % Fragmentation (at any interior point of linear chain)
    });

%% Add reactions

stoichTempl = zeros(size(Model.species, 1), 1);

% Primary nucleation reaction
stoichPN = stoichTempl;
stoichPN(1) = -1 * nc; % n_c monomers are consumed to make the nucleus
stoichPN(nc) = 1; % One oligomer of size n_c is produced
Model = Model.addReaction(['kn*x1^', num2str(nc)], stoichPN);

% Fragmentation reactions
% Any nucleus of length [nc, MAX_NUCLEUS_SIZE] can be fragmented at any
% interior point in the chain. This will split the existing oligomer into
% two oligomers whose combined length will equal that of the original 
% chain.

for l = nc : MAX_NUCLEUS_SIZE
    for i = 1 : floor(l/2)
        stoichF = stoichTempl;
        stoichF(i) = 1; % One oligomer of length i is produced
        stoichF(l-i) = 1; % One oligomer of length l-i is produced
        stoichF(l) = -1; % One oligomer of length l is consumed

        % The propensity will depend on the number of locations in the
        % chain where an equivalent fragmentation could occur. For an
        % exact bifurcation, there is only a single suitable point, e.g.
        % only one cut will split a 6-mer into two 3-mers. For inexact
        % bifurcations, e.g. a 7-mer into one 2- and one 5-mer, there are
        % two cuts:

        if i < (l/2)
            propensity = ['2*kminus*x', num2str(l)];
        else
            propensity = ['kminus*x', num2str(l)];
        end

        Model = Model.addReaction(propensity, stoichF);
    end % for i = 1 : floor(l/2) ...
end % for l = nc : MAX_NUCLEUS_SIZE ...

% Secondary nucleation reactions
% Any nucleus of length [nc, MAX_NUCLEUS_SIZE - n2] can grow via the 
% attachment of another nucleus at any interior point in the chain.
% As in the case of fragmentation reactions, the number of possible
% attachment sites depends on the length of the existing chain. 

for l = nc : (MAX_NUCLEUS_SIZE - n2)
    stoichSN = stoichTempl;      
    stoichSN(l) = -1; % One oligomer of length l is consumed
    stoichSN(1) = -n2; % n2 monomers are consumed
    stoichSN(l+n2) = 1; % One oligomer of length l+n2 is produced

    % The propensity will depend on the number of locations in the
    % chain where an equivalent attachment could occur. We assume that
    % chains are linear, before and after; therefore, there is
    % symmetry throughout, and the number of equivalent attachment 
    % sites equals the number of interior monomers, e.g. a 7-mer will
    % have five interior monomers.

    propensity = ['k2*(x1^', num2str(n2), ')*x', num2str(l), ...
        '*', num2str(l-2)];

    Model = Model.addReaction(propensity, stoichSN);
end % for l = nc : (MAX_NUCLEUS_SIZE - n2) ...

% Elongation reactions
% Any nucleus of length [nc, MAX_NUCLEUS_SIZE) can grow via the 
% attachment of a monomer at either endpoint of the chain. Regardless of
% the length of the existing chain, there are exactly two possible points
% of attachment (since the existing nucleus is at minimum a dimer, it has
% two distinct ends).

for l = nc : (MAX_NUCLEUS_SIZE - 1)
    stoichE = stoichTempl;      
    stoichE(l) = -1; % One oligomer of length l is consumed
    stoichE(1) = -1; % One monomer is consumed
    stoichE(l+1) = 1; % One oligomer of length l+1 is produced

    propensity = ['2*kplus*x1*x', num2str(l)];

    Model = Model.addReaction(propensity, stoichE);
end % for l = nc : (MAX_NUCLEUS_SIZE - 1) ...

Model.summarizeModel

%% Generate Equations for Propensities
Model = Model.formPropensitiesGeneral('FibrilFormation');

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
JUncertain = find(diag(Vals) < 1e-9);
figure
bar(Vecs(:,JUncertain))