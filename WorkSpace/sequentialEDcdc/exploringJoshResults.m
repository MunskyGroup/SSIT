firstExperiment = zeros(1,12);
firstExperiment([1,6,12]) = round(120/3);

load IterativeExperimentResults_DUSP1_real_FIMopt_4_10rd
FIM_fim = FIMcurrentExptTrueSaved{end};

expFIM = firstExperiment;
for i = 1:9
    expFIM = expFIM + exptDesigns{i};
end

load IterativeExperimentResults_DUSP1_real_random_4_10rd
FIM_rand = FIMcurrentExptTrueSaved{end};
expRAND = firstExperiment;
for i = 1:9
    expRAND = expRAND + exptDesigns{i};
end

bar([expFIM',expRAND'])


%%
ModelTrue = SSIT; % Rename SSIT code to make changes with out changing the original code
ModelTrue.species = {'x1';'x2'}; % x1: number of on alleles,  x2: mRNA
ModelTrue.initialCondition = [0;0]; % Set initial conditions
ModelTrue.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'}; % Set the propensity functions of the 4 reactions
% Associated stoichiometry of the four reactions
stoich = [1,0;   % Gene allele switching to on state
    -1,0;   % Gene allele switching to off state
    0,1;   % mRNA production
    0,-1]; % mRNA degredation
ModelTrue.stoichiometry = transpose(stoich);
% Input expression for time varying signals
ModelTrue.inputExpressions = {'IGR','kcn0/knc+(t>=0)*kcn1/(r1-knc)*(exp(-knc*t)-exp(-r1*t))'};
% Defining the values of each parameter used

ModelTrue.parameters=load('SGRS_model_v1.mat').SGRS_Model.parameters; % True model parameters after fitting to data

ModelTrue.fspOptions.initApproxSS = true;
        
fitParameters = [1:4];
timeDUSP1 = [0 10 20 30 40 50 60 75 90 120 150 180];
nT = length(timeDUSP1);
ModelTrue.tSpan = timeDUSP1;
ModelTrue.fspOptions.bounds(3:4) = [2.1,100];
fimMetric = 'TR1:4';

ModelTrue.pdoOptions.unobservedSpecies = {'x1'};

muLog10Prior = [-1 -1 -1 -2 -2 -2 -2 -2];
sigLog10Prior = [2 2 2 2 2 2 2 2];
ModelTrue.fittingOptions.modelVarsToFit = fitParameters;

% muLog10Prior = muLog10Prior(ModelTrue.fittingOptions.modelVarsToFit);
% sigLog10Prior = sigLog10Prior(ModelTrue.fittingOptions.modelVarsToFit);
% ModelTrue.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);

nSamplesMH = 5000; % Number of MH Samples to run

intuitiveDesign = [20 0 20 0 20 0 20 0 20 0 0 20;
    30 0 0 30 0 0 30 0 0 0 30 0;
    0 0 0 40 0 0 40 0 0 40 0 0;
    0 0 0 40 0 0 0 40 0 40 0 0;
    0 0 0 0 40 0 40 0 0 40 0 0;
    0 0 40 0 40 0 0 40 0 0 0 0;
    40 0 0 0 40 0 0 0 0 0 0 0;
    0 0 0 40 0 0 0 40 0 40 0 0;
    0 0 0 0 40 0 40 0 0 40 0 0;
    80 0 0 0 40 0 0 0 0 0 0 0];

% Total number of new cells allowed in each experiment
numCellsPerExperiment = 120;
ModelTrue = ModelTrue.formPropensitiesGeneral('DUSP1_tmp',true);
ModelTrue.fspOptions.fspTol = 1e-6;
ModelTrue.fspOptions.verbose = true;
[ModelSolution,ModelTrue.fspOptions.bounds] = ModelTrue.solve;

%%
FIM = ModelTrue.computeFIM([],'log');
FIMTotal_FIMExpt = ModelTrue.totalFim(FIM,expFIM);
FIMTotal_RandExpt = ModelTrue.totalFim(FIM,expRAND);

det_FimExpt = det((FIMTotal_FIMExpt{1}(1:4,1:4))^-1)
det_RandExpt = det((FIMTotal_RandExpt{1}(1:4,1:4))^-1)

% maxEV_fim = max(eig((FIMTotal_FIMExpt{1}(1:4,1:4))^-1))
% maxEV_rand = max(eig((FIMTotal_RandExpt{1}(1:4,1:4))^-1))

detFimTrueInv_fim_Josh = det((FIM_fim{1}(1:4,1:4))^-1)
detFimTrueInv_rand_Josh = det((FIM_rand{1}(1:4,1:4))^-1)

optimum = ModelTrue.optimizeCellCounts(FIM,sum(expRAND),'TR1:4');
FIMTotal_Optimum = ModelTrue.totalFim(FIM,optimum);
det_Fim_Opt = det((FIMTotal_Optimum{1}(1:4,1:4))^-1)




