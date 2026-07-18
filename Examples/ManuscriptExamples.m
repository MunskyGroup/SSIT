% Install all SSIT tools and packages
% install

% Load STL1_4state model (see SI for manual creation)
STL1_4state = SSIT('STL1_4state');
% Solve ODEs:
STL1_4state = STL1_4state.solve(solver='ODE');
% STL1_4state.solutionScheme = 'ODE';
% STL1_4state.Solutions = STL1_4state.solve; 

% Plot ODE solutions for mRNA:
STL1_4state.plotODE(speciesNames={'mRNA'}, timeVec=STL1_4state.tSpan, ...
    lineProps={'linewidth',4}, TitleFontSize=26, Title='4-state STL1 (mRNA)', ...
    AxisLabelSize=20, TickLabelSize=20, LegendFontSize=20, LegendLocation='east',...
    Colors=[0.23,0.67,0.20], XLabel='Time', YLabel='Molecule Count')
% Plot ODE solutions for the four gene states:
STL1_4state.plotODE(speciesNames=STL1_4state.species(1:4), ...
    Title='4-state STL1 (gene states)', YLabel='Molecule Count')

%%

% Set the number of simulations performed per experiment:
STL1_4state.ssaOptions.nSims=100;
% Use a negative initial time to allow model to equilibrate to an initial
% steady state prior to starting the subsequent simulation (burn-in):
STL1_4state.tSpan = [-100,STL1_4state.tSpan];
% Set the initial time:
STL1_4state.initialTime = -100;
% Run SSA:
STL1_4state = STL1_4state.solve(solver='SSA');

%%
% Set FSP 1-norm error tolerance:
STL1_4state.fspOptions.fspTol = 1e-4;

% Guess initial bounds on FSP StateSpace:
STL1_4state.fspOptions.bounds = [0;0;0;0;0;1;1;1;1;200];

% Approximate the steady state for the initial distribution:
STL1_4state.fspOptions.initApproxSS = true;
STL1_4state.tSpan = [0:5:60];
STL1_4state.initialTime = 0;

% Solve Model:
STL1_4state = STL1_4state.solve(solver='FSP');

%%
% Solve for time for mRNA to reach 100:
STL1_4state.fspOptions.escapeSinks.f = {'mRNA'};
STL1_4state.fspOptions.escapeSinks.b = 100;
STL1_4state = STL1_4state.solve(solver='FSP');

%%
% Solve the sensitivity problem:
STL1_4state = STL1_4state.solve(solver='fspSens');

%%

% Specify time points and numbers of cells for experiments
STL1_4state.tSpan = [0,1,2,4,6,8,10:5:55];
cellCounts = 1000*ones(1,16);

% Compute FIMs using FSP sensitivity results (for all parameters)
fims_full = STL1_4state.computeFIM(scale='log',observed='mRNA'); % full FIM

% Compute FIMs using FSP sensitivity results (Hog1 parameters fixed)
fims_free = STL1_4state.computeFIM(scale='log',observed='mRNA',...
    freePars=freePars); % FIM sub matrix

% Evaluate the provided experiment design (in "cellCounts")
% to get the FIM:
[fimTotal_full,mleCovEstimate_full,fimMetrics_full] = ...
STL1_4state.evaluateExperiment(fims_full,cellCounts)

[fimTotal_free,mleCovEstimate_free,fimMetrics_free] = ...
STL1_4state.evaluateExperiment(fims_free,cellCounts)

%%

% Experiment Design: Find the FIM-based designs for total of 16k cells
% Compute the optimal number of cells from the FIM results using different
% design criteria:
nTotal = 16000;
nC_Dcov = STL1_4state.optimizeCellCounts(fims_full,nTotal,'D-cov');
nC_Trace = STL1_4state.optimizeCellCounts(fims_full,nTotal,'Trace');
nC_Dopt8 = STL1_4state.optimizeCellCounts(fims_full,nTotal,'D-opt-sub[1:8]');
nC_Dopt13 = STL1_4state.optimizeCellCounts(fims_full,nTotal,'D-opt-sub[9:13]');
nC_Dopt18 = STL1_4state.optimizeCellCounts(fims_full,nTotal,'D-opt-sub[14:18]');

%%
STL1_4state = STL1_4state.loadData('data/filtered_data_2M_NaCl_Step.csv',...
{'mRNA','RNA_STL1_total_TS3Full'},...
{'Replica',1;'Condition','0.2M_NaCl_Step'});

%%
% An example of fitting model parameters to loaded data using maximum likelihood
% estimation (MLE) is shown below.

% Set fitOptions, with the maximum allowable number of iterations to fit:
fitOptions = optimset('Display','iter','MaxIter',2000);

% Specify parameters to fit (all but Hog1 input parameters, which are fixed):
STL1_4state.fittingOptions.modelVarsToFit = 1:13;

% Search to Find the MLE:
STL1_4state = STL1_4state.maximizeLikelihood(fitOptions=fitOptions);

%%
% Specify Prior as log-normal distribution with wide uncertainty.
mu_log10 = [0.5,2,5,3.5,-0.4,1,0.2,0.4,-0.5,-1.3,-0.1,2,0.5]; % Prior log-mean
sig_log10 = 2*ones(1,13); % Prior log-standard deviation

% Log Prior (neglecting the normalization constant):
STL1_4state.fittingOptions.logPrior = ...
@(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Set MH proposal options, adjust aiming for acceptance ratio around 0.3-0.4:
proposalWidthScale = 0.01;
MHOptions.proposalDistribution = @(x)x+proposalWidthScale*randn(size(x));

% Set MH runtime options (number of samples, burnin, thinning, etc.):
MHOptions.numberOfSamples = 20; MHOptions.burnin = 5; MHOptions.thin = 2;
% Run Metropolis-Hastings:
STL1_4state = STL1_4state.maximizeLikelihood(fitOptions=MHOptions, fitAlgorithm='MetropolisHastings');

%%
% Select models to include in multimodel:
ModelNames = {'Model_DUSP1','Model_RUNX1','Model_BIRC3'};
for i = 1:3
    % Load and prepare codes for each model
    Models{i} = SSIT(['seqModels/',ModelNames{i}]);
    % Set free parameters
    Models{i}.fittingOptions.modelVarsToFit = 1:9;
    % Assign shared (1,2) and unshared (3,..,9) parameters:
    ParInds{i} = [1,2,(3:9)+7*(i-1)];
end
Models{1}.formPropensitiesGeneral;
% Set a constraint on model parameters:
Constraint = @(x) -var(log10([x(3:9);x(7*1+(3:9));x(7*2+(3:9))]));
% Create and initialize multimodel:
combinedModel = SSITMultiModel(Models, ParInds, Constraint);
combinedModel = combinedModel.initializeStateSpaces();
% Fit the multimodel:
fitOptions = optimset('Display','iter','MaxIter',1000);
combinedModel = combinedModel.maximizeLikelihood(fitOptions=fitOptions);

%% Call Pipeline to Fit scRNA-seq Models
% Specify pipeline to apply to model and arguments
% ("../../SSIT/src/exampleData/examplePipelines/fittingPipelineExample.m")
Pipeline = 'fittingPipelineExample';
pipelineArgs.maxIter = 1000; pipelineArgs.display = 'iter';
pipelineArgs.makePlot = false; pipelineArgs.nRounds = 1;
% Launch cluster jobs for all genes.
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
TAB = readtable(DataFileName);
geneNames = fields(TAB);
for iGene = 2:length(geneNames)-4
    modelName = ['Model_',geneNames{iGene}];
    saveName = ['seqModels/',modelName];
    logfile = ['logFiles/log',modelName];
    if iGene==2
        load(saveName)
        eval([modelName,'.formPropensitiesGeneral;']);
    end
    cmd = SSIT.generateCommandLinePipeline(saveName,modelName,Pipeline,...
        pipelineArgs=pipelineArgs,saveFileOut=saveName, ...
        logFile=logfile,runNow=true,runOnCluster=true);
end
