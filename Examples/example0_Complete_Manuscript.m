clear

%% 2.1
% Add SSIT source codes to Matlab search path.
addpath(genpath('../src'));
STL1_4state = SSIT('Empty');
STL1_4state.species = {'g1'; 'g2'; 'g3'; 'g4'; 'mRNA'};
STL1_4state.initialCondition = [1;0;0;0;0];

%% 2.2
STL1_4state.inputExpressions = ...
    {'Hog1',['A*(((1-(exp(-r1*(t-t0))))*',...
    'exp(-r2*(t-t0)))/(1+((1-(exp(-r1*(t-t0))))*',...
    'exp(-r2*(t-t0)))/M))^n*(t>t0)']};

%% 2.3
STL1_4state.stoichiometry = [-1, 1, 0, 0, 0, 0, 0;... % gene state 1
1,-1,-1, 1, 0, 0, 0;... % gene state 2
0, 0, 1,-1,-1, 1, 0;... % gene state 3
0, 0, 0, 0, 1,-1, 0;... % gene state 4
0, 0, 0, 0, 0, 0, 1]; 

newReaction.propensity = 'dr*mRNA';
newReaction.stoichiometry = {'mRNA',-1};
newReaction.parameters = {'dr',1};
STL1_4state = STL1_4state.addReaction(newReaction);

%% 2.4
STL1_4state.propensityFunctions = {...
    'k12*g1'; '(max(0,k21o*(1-k21i*Hog1)))*g2';...
    'k23*g2'; 'k32*g3';...
    'k34*g3'; 'k43*g4';...
    'kr*g4' ; 'dr*mRNA'};

%% 2.5
STL1_4state.parameters = ({'t0',5.8; 'k12',90;
    'k21o',1e+03; 'k21i',1;
    'k23',5e+02; 'k34',5;
    'k32',1000; 'k43',200;
    'dr',1; 'kr',2500;
    'r1',35; 'r2',0.5;
    'A',3; 'M',25;
    'n',0.1});

%% 2.6
STL1_4state.tSpan = linspace(0,50,101);

%% 2.7
STL1_4state.summarizeModel

%% 2.8
save('example_1_CreateSSITModels','STL1_4state')
load('example_1_CreateSSITModels.mat')

%% 2.9
% Set solution scheme to 'ODE':
STL1_4state.solutionScheme = 'ODE';

% Compile and store the given reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_ODE');

% Solve ODEs:
[~,~,STL1_4state] = STL1_4state.solve;

% Plot ODE solutions for mRNA:
STL1_4state.plotODE(STL1_4state.species(5), STL1_4state.tSpan,...
    {'linewidth',4}, Title='4-state STL1 (mRNA)', TitleFontSize=24,...
    LegendLocation='east',Colors=[0.23,0.67,0.20], YLabel='Molecule Count')

% Plot ODE solutions for the 4 gene states:
STL1_4state.plotODE(STL1_4state.species(1:4), STL1_4state.tSpan,...
    {'linewidth',4}, Title='4-state STL1 (gene states)', YLim=[0,1],...
    TitleFontSize=24, LegendLocation='east', YLabel='Molecule Count')

%% 2.10
% Compile and store the given reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_moments');

% Set solution scheme and solve.
STL1_4state.solutionScheme = 'moments';
[~,~,STL1_4state] = STL1_4state.solve;

STL1_4state.plotMoments(STL1_4state.Solutions, STL1_4state.species(5),...
    "meansanddevs", STL1_4state.tSpan, [], {'linewidth',4},...
    Title='4-state STL1 (mRNA)', TitleFontSize=24,...
    LegendLocation='northeast', Colors=[0.23,0.67,0.20],...
    YLabel='Molecule Count')

STL1_4state.plotMoments(STL1_4state.Solutions, STL1_4state.species(1:4),...
    "meansanddevs", STL1_4state.tSpan, [], {'linewidth',4},...
    Title='4-state STL1 (gene states)', TitleFontSize=24,...
    LegendLocation='northeast', YLabel='Molecule Count')

%% 2.11
% Set solution scheme to SSA:
STL1_4state.solutionScheme = 'SSA';

% Set the number of simulations performed per experiment
STL1_4state.ssaOptions.Nsims = 100;

% Use a negative initial time to allow model to equilibrate (burn-in)
% to an initial steady state prior to starting the subsequent simulation:
STL1_4state.tSpan = [-100,STL1_4state.tSpan];

% Set the initial time:
STL1_4state.initialTime = STL1_4state.tSpan(1);

% Run iterations in parallel with multiple cores, or execute serially:
STL1_4state.ssaOptions.useParallel = true;

% Compile and store the given reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_SSA');

% Run SSA:
[~,~,STL1_4state] = STL1_4state.solve;

% Plot SSA trajectories and means (mRNA):
STL1_4state.plotSSA('all',100,STL1_4state.species(5),{'linewidth',4}, ...
    HistTime=20, Title="4-state STL1 (mRNA)", MeanOnly=false,...
    YLabel='Molecule Count', Colors=[0.23,0.67,0.20], TitleFontSize=24,...
    LegendLocation='northeast');

%% FSP
% Ensure the solution scheme is set to FSP (default):
STL1_4state.solutionScheme = 'FSP';

% Set FSP 1-norm error tolerance:
STL1_4state.fspOptions.fspTol = 1e-4;

% Guess initial bounds on FSP StateSpace:
STL1_4state.fspOptions.bounds = [0,0,0,0,0,1,1,1,1,200];

% Approximate the steady state for the initial distribution:
STL1_4state.fspOptions.initApproxSS = true;
STL1_4state.tSpan = [0:5:60];
STL1_4state.initialTime = 0;

% Compile and store the given reaction propensities: 
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_FSP');

% Solve Model:
[~,~,STL1_4state] = STL1_4state.solve;

% Plot means and standard deviations:
STL1_4state.plotFSP(STL1_4state.Solutions, STL1_4state.species(5),...
    'meansAndDevs', [], [], {'linewidth',4}, YLabel='Molecule Count',...
    Title='4-state STL1 (mRNA)', TitleFontSize=24, LegendFontSize=18,...
    LegendLocation='northeast', Colors=[0.23,0.67,0.2]);

% Plot marginal distributions:
STL1_4state.plotFSP(STL1_4state.Solutions, STL1_4state.species(5),...
    'marginals', [1,3,7,11], [], {'linewidth',3},...
    Colors=[0.23,0.67,0.2], XLim=[0,100])

%% FSP Escape Times
% Make copy of original model
STL1_4state_escape = STL1_4state;

% Set the initial populations:
STL1_4state_escape.initialCondition = [1;0;0;0;0];

% Turn off steady state initial condition
STL1_4state_escape.fspOptions.initApproxSS = false;

% Set the times at which distributions will be computed:
STL1_4state_escape.tSpan = linspace(0,60,100);
STL1_4state_escape.initialTime = 0;

% Solve for time to reach mRNA=100:
STL1_4state_escape.fspOptions.escapeSinks.f = {'mRNA'};
STL1_4state_escape.fspOptions.escapeSinks.b = 100;
[~,~,STL1_4state_escape] = STL1_4state_escape.solve;

% Plot the CDF and PDF
STL1_4state_escape.plotFSP(STL1_4state_escape.Solutions, [], "escapeTimes",...
    [], [], {'linewidth',3}, TitleFontSize=24, Title="4-state STL1 (mRNA)",...
    Colors=[0.23,0.67,0.2], LegendLocation="southeast", XLim=[0,50]);

%% Solve FSP sensitivities
% Set solution scheme to FSP sensitivity:
STL1_4state.solutionScheme = 'fspSens';

% Compile and store the given reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_sens');

% Solve the sensitivity problem:
[~,~,STL1_4state] = STL1_4state.solve;

% Choose time at which to make sensitivity plots
index_PlotTime = 6;  

% Plot the results from the sensitivity analysis
STL1_4state.plotFSP(STL1_4state.Solutions,...
    STL1_4state.species(5), 'sens', index_PlotTime, [], {'linewidth',3}, ...
    Colors=[0.23,0.67,0.2], AxisLabelSize=15, TickLabelSize=12, ...
    XLim=[0,10], Title="4-state STL1 (t=25)", TitleFontSize=24)

%% FIM Analysis

% Set unobservable species:
STL1_4state.pdoOptions.unobservedSpecies = '1:4';

% Add lognormal prior (so FIM is invertible). Prior log-standard deviation:
sig_log10 = 2*ones(1,15);

% Compile and store the given reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_FIM');

% Compute FIMs using FSP sensitivity results.
fimResults = STL1_4state.computeFIM;

% Specify how many cells are to be measured at each time:
cellCounts = 1000*ones(size(STL1_4state.tSpan));

% Evaluate the provided experiment design (in "cellCounts")
% and produce an array of FIMs (one for each parameter set):
[STL1_4stateTotalFIM, STL1_4state_mleCovEstimate, STL1_4stateMetrics] = ...
STL1_4state.evaluateExperiment(fimResults, cellCounts, diag(sig_log10.^2));

% Plot the FIMs:
STL1_4state.plotFIMResults(STL1_4stateTotalFIM, STL1_4state.parameters,...
    [STL1_4state.parameters{:,2}], EllipsePairs=[1 6; 9 10; 8 10; 8 9]);

%% Experiment Design

% Find the FIM-based designs for a total of 1000 cells 

% Compile and store propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_design');

% Compute the optimal number of cells from the FIM results using different 
% design criteria:  `Trace' maximizes the trace of the FIM; 
% `DetCovariance' minimizes the expected determinant of MLE covariance; 
% `Smallest Eigenvalue' maximizes the smallest e.val of the FIM; and 
% `TR[$<i_1>,<i_2>$,...]' maximizes the determinant of the FIM for the 
% specified indices.  The latter is shown for different parameter 
% combinations, where `Tr[9:10]' are the mRNA-specific parameters `dr' and 
% `kr' (degradation and transcription, respectively).  All other parameters 
% are assumed to be known and fixed.
nCol = sum(cellCounts);
nTotal = nCol(1);
nCellsOpt_detCov = ...
    STL1_4state.optimizeCellCounts(fimResults,nTotal,'DetCovariance');
nCellsOpt_trace = ...
    STL1_4state.optimizeCellCounts(fimResults,nTotal,'Trace');
nCellsOpt_tr = ...
    STL1_4state.optimizeCellCounts(fimResults,nTotal,'tr[1:10]');
nCellsOpt_tr1 = ...
    STL1_4state.optimizeCellCounts(fimResults,nTotal,'tr[11:15]');
nCellsOpt_trR = ...
    STL1_4state.optimizeCellCounts(fimResults,nTotal,'tr[9:10]');

%% Make a bar chart to compare the different designs
% Find which x positions correspond to time=30 and time=60 for off-setting:
x = 1:13;
t = STL1_4state.tSpan;               
idx = ismember(t, [30 60]);        

% Build custom x-locations for series that need separation
x_m02 = x;  x_m03(idx) = x_m02(idx) - 0.5;
x_p02 = x;  x_p02(idx) = x_p02(idx) + 0.5;

f = figure;
bar(x,      nCellsOpt_trace,  0.5); hold on
bar(x_m02,  nCellsOpt_trR,    0.5);
bar(x_p02,  nCellsOpt_tr1,    0.5);
bar(x,      nCellsOpt_detCov, 0.5);
bar(x_m02,  nCellsOpt_tr,     0.5);

set(gca,'XTick',x,'XTickLabel',t,'FontSize',16)
title('4-state STL1 (FIM Optimal Designs)','FontSize',24)
xlabel('Time (min)','FontSize',20)
ylabel('Number of cells','FontSize',20)
legend('Trace,Tr[1] Designs','Tr[9:10] Design','Tr[11:15] Design', ...
       'DetCov,\lambda Designs','Tr[1:10] Design','Location','northeast')


fimOpt_tr = STL1_4state.evaluateExperiment(fimResults,nCellsOpt_tr,...
                                           diag(sig_log10.^2));

% Plot the FIMs for 'tr[1:10]':
STL1_4state.plotFIMResults(fimOpt_tr, STL1_4state.parameters,...
    [STL1_4state.parameters{:,2}], EllipsePairs=[1 6; 9 10; 8 10; 8 9]);


%% Load and Plot Data
% Note: Ensure the search path is correct on the local machine: 
STL1_4state = ...
    STL1_4state.loadData('data/filtered_data_2M_NaCl_Step.csv',...
    {'mRNA','RNA_STL1_total_TS3Full'},...
    {'Replica',1;'Condition','0.2M_NaCl_Step'});

% This plot is unnecessary, as the model parameters have not been fit to
% the data yet. However, it illustrates the improvement to come later:
STL1_4state.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);


%% Find MLE
% Maximum allowable number of iterations to fit, etc.:
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them):
STL1_4state.fittingOptions.modelVarsToFit = [1:15];

% Search to Find the MLE:
[~,~,~,STL1_4state] = STL1_4state.maximizeLikelihood([],fitOptions);

% Make plots of the parameter fits from the MLE:
STL1_4state.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

% Note: Should see an MLE of -26015 at the end:
% Exiting: Maximum number of iterations has been exceeded
%          - increase MaxIter option.
%          Current function value: 26014.985315

%% Specify Bayesian Prior and fit

% Make a copy of our 4-state STL1 model:
STL1_4state_MH = STL1_4state;

% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [0.8,0,0.3,1.2,-1,1,3.5,0,2,3,0.5,3.5,3,-1,4];

% Prior log-standard deviation:
sig_log10 = 2*ones(1,15);  

% Prior:
STL1_4state_MH.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose parameters to search:
STL1_4state_MH.fittingOptions.modelVarsToFit = [1:15]; 

% Create first parameter guess:
STL1_4state_MH_pars = [STL1_4state_MH.parameters{:,2}];      

% Iterating between MLE and MH
%  Running a few rounds of MLE and MH together may improve convergence:

STL1_4state_MH.parameters(:,2) = num2cell(STL1_4state_MH_pars);
for i=1:2
    % Maximize likelihood:
    STL1_4state_MH_pars = STL1_4state_MH.maximizeLikelihood([]);    
    % Update parameters in the model:
    STL1_4state_MH.parameters(:,2) = num2cell(STL1_4state_MH_pars);

    % Run Metropolis-Hastings    
    proposalWidthScale = 0.005;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

    % Set MH runtime options (number of samples, burnin, thin, etc.):
    MHOptions.numberOfSamples = 2000;
    MHOptions.burnin = 200;
    MHOptions.thin = 2;

    % Run Metropolis-Hastings: 
    [STL1_4state_MH_pars,~,STL1_4state_MH_MHResults] = ...
        STL1_4state_MH.maximizeLikelihood([], MHOptions,...
        'MetropolisHastings');
    
    % Store MH parameters in model:
    STL1_4state_MH.parameters([1:15],2) = ...
        num2cell(STL1_4state_MH_pars);
end
% STL1_4state_MH.plotMHResults(STL1_4state_MH_MHResults,...
%     paramSelect={'t0','k12','k21o','k21i','kr','dr'});
STL1_4state_MH.plotMHResults(STL1_4state_MH_MHResults);

STL1_4state_MH.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count', ProbXLim = [0 80],...
    LegendLocation='northeast', LegendFontSize=12,...
    TimePoints=[0 8 10 15 30 55]);

% TODO - Need a wrapper for the MH Plotting functions to allow easier
% choice of which scatter plots to show.

%% ABC
STL1_4state_ABC = STL1_4state;
% Set up a prior over parameters (logPriorLoss)
logPriorLoss = @(x)sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose loss function for ABC (default: 'cdf_one_norm'):
lossFunction = 'cdf_one_norm';

% Set ABC / MCMC options
ABCoptions = struct('numberOfSamples',5,'burnIn',0,'thin',1,...
    'proposalDistribution',@(x)x+0.01*randn(size(x)));

% Compile and store reaction propensities:
STL1_4state_ABC = STL1_4state_ABC.formPropensitiesGeneral('STL1_4state_ABC');

% Run ABC search
[~, ~, ~, STL1_4state_ABC] = STL1_4state_ABC.runABCsearch([],...
    lossFunction, logPriorLoss, ABCoptions);

STL1_4state_ABC.plotABC(STL1_4state_ABC.Solutions.ABC);
%%
if isfield(ResultsABC, 'mhSamples')
    parChain = ResultsABC.mhSamples;   % size: [numberOfSamples x nPars] 
    nPars    = size(parChain, 2);

    figure;
    for k = 1:nPars
        subplot(ceil(nPars/2), 2, k);
        histogram(parChain(:,k), 40, 'Normalization', 'pdf');
        hold on;
        xline(parsABC(k), 'r', 'LineWidth', 1.5);
        xline(STL1_4state_MH_pars(k), 'b', 'LineWidth', 1.5);
        title(sprintf('Parameter %d', k));
        xlabel('\theta_k');
        ylabel('Posterior density (approx.)');
    end
    sgtitle('ABC posterior marginals (approximate)');
else
    warning('ResultsABC.mhSamples not found.');
end


%% Cross Validation
% Specify datafile name and species linking rules:
DataFileName = 'data/filtered_data_2M_NaCl_Step.csv';
LinkedSpecies = {'mRNA','RNA_STL1_total_TS3Full'};

% Set the global conditions (e.g., fit data at times before 50 min.):
ConditionsGlobal = {[],[],'TAB.time<=50'};

% Split up the replicas to be separate:
ConditionsReplicas = {'TAB.Replica==1';'TAB.Replica==2'};

% Specify constraints on rep-to-rep parameter variations. Here, we specify
% that there is an expected 0.1 log10 deviation expected in some parameters
% and smaller in others. No deviation at all is indicated by 0.
Log10Constraints = ...
    [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.02,0.02,0.02,0.02,0.02,0.1,0.1];

% Compile and store reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state');

% Create full model:
CrossValidationModel = SSITMultiModel.createCrossValMultiModel(...
    STL1_4state, DataFileName, LinkedSpecies, ConditionsGlobal,...
    ConditionsReplicas, Log10Constraints);

% Run the model fitting routines:
[~,~,~,CrossValidationModel] = ...
    CrossValidationModel.maximizeLikelihood([],fitOptions);

% Make a figure to explore how much the parameters changed between replicas:
fignum = 12; useRelative = true;
CrossValidationModel.compareParameters(fignum,useRelative);


%% Multi-Model
% This looks very similar to the previous cross-validation model. I recommend
% combining it with the above.

%% Model Reduction
% None of the current model reductions are meant for use in time varying
% problems, so I doubt that they would work for the Hog Model.  Also, with
% my recent changes, the Hog1 model is a lot faster than before.
% I recommend removing this section and demonstrate on a different model in
% the SI.

% Make a copy of the STL1 model to set up for model reduction:
STL1_ModRed = STL1_4state;
STL1_ModRed.fspOptions.initApproxSS = true;
STL1_ModRed.modelReductionOptions.useModReduction = true;
STL1_ModRed.fspOptions.fspTol = inf;
STL1_ModRed.modelReductionOptions.reductionType = 'Logarithmic State Lumping';
STL1_ModRed.modelReductionOptions.reductionOrder = 30;
[STL1_ModRed,fspSets] = STL1_ModRed.computeModelReductionTransformMatrices;

tic
redSoln = STL1_ModRed.solve(fspSets.stateSpace);
redModelSolveTime = toc

tic
fullSoln = STL1_4state.solve(fspSets.stateSpace);
fullModelSolveTime = toc

% Plot the full and reduced FSP solutions:
STL1_ModRed.plotFSP(redSoln,...
    STL1_ModRed.species(5), 'meansAndDevs', [], [], {'linewidth',4},...
    Title='4-state STL1 (FSP Reduced)', TitleFontSize=24,...
    XLabel='Time', Colors=[0.23,0.67,0.2], YLabel='Molecule Count',...
    LegendFontSize=15, LegendLocation='northeast',YLim=[0,40]);

% Plot the full and reduced FSP solutions:
STL1_4state.plotFSP(fullSoln,...
    STL1_ModRed.species(5), 'meansAndDevs', [], [], {'linewidth',4},...
    Title='4-state STL1 (FSP Full)', TitleFontSize=24,...
    XLabel='Time', Colors=[0.23,0.67,0.2], YLabel='Molecule Count',...
    LegendFontSize=15, LegendLocation='northeast',YLim=[0,40]);

%% Hybrid Models
% This extended model is based on that published here: 
% https://www.cell.com/fulltext/S0092-8674(09)00508-X?large_figure=true
% The parameters are taken from the paper, except for k_h, which is set to
% match the approximate magnitude of Hog1 in the above model.

% Make copy of previous model
STL1_4state_Extended = STL1_4state;

% Remove the input Hog1(t) to be replaced upstream reactions
STL1_4state_Extended.inputExpressions = {};

% Add reaction species, stoichiometries, propensities, and parameters. When
% new species are detected, their initial conditions are automatically set
% to zero.
STL1_4state_Extended = STL1_4state_Extended.addReaction(struct(...
    'propensity','k_h*NaCl*(t>t0)/(1+(G/NaCl))','stoichiometry',{{'Hog1',1}},...
    'parameters',{{'k_h',10;'NaCl',0.2}}));
STL1_4state_Extended = STL1_4state_Extended.addReaction(struct(...
    'propensity','gamma_h*Hog1','stoichiometry',{{'Hog1',-1}},...
    'parameters',{{'gamma_h',0.369}}));
STL1_4state_Extended = STL1_4state_Extended.addReaction(struct(...
    'propensity','alpha_d*Hog1','stoichiometry',{{'D',1}},...
    'parameters',{{'alpha_d',0.0106}}));
STL1_4state_Extended = STL1_4state_Extended.addReaction(struct(...
    'propensity','D+alpha_i*NaCl*(t>t0)/(1+(G/NaCl))','stoichiometry',{{'G',1}},...
    'parameters',{{'alpha_i',0.0806}}));
STL1_4state_Extended = STL1_4state_Extended.addReaction(struct(...
    'propensity','gamma_g*G','stoichiometry',{{'G',-1}},...
    'parameters',{{'gamma_g',0.119}}));

% Set model to emply hybrid solving approach:
STL1_4state_Extended.useHybrid = true;

% Define which species will be replaced with upstream ODEs:
STL1_4state_Extended.hybridOptions.upstreamODEs = {'Hog1','D','G'};

% Compile and store reaction propensities:
STL1_4state_Extended = ...
    STL1_4state_Extended.formPropensitiesGeneral('STL1_4state_Extended');

% Solve hybrid ODE-FSP model
[~,~,STL1_4state_Extended] = STL1_4state_Extended.solve;

% Plot the results
STL1_4state_Extended.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

% TODO - need to update the plotting function to also allow for plots of
% the upstream ODEs. This might be in the codes already.

%%
% close all
% % Plot ODE solutions for mRNA:
% STL1_4state_Extended.plotODE(STL1_4state_Extended.species(6),...
%     STL1_4state_Extended.tSpan, {'linewidth',4},...
%     Title='4-state STL1 (mRNA)', TitleFontSize=24,...
%     AxisLabelSize=18, TickLabelSize=18, LegendFontSize=15,...
%     LegendLocation='east', Colors=[0.23,0.67,0.20],...
%     XLabel='Time', YLabel='Molecule Count')


% The following can be used to generate the plot of the hog signal
% from the previous model
% Compare to the other function for the hog model.
% hold on
% t = linspace(0,60,100);
% % Need to adjust these parameters to match the actual data.
% r1=STL1_4state.parameters{11,2};
% r2=STL1_4state.parameters{12,2};
% A=STL1_4state.parameters{13,2}; 
% M=STL1_4state.parameters{14,2};
% n=STL1_4state.parameters{15,2};
% t0=STL1_4state.parameters{1,2};
% 
% Hog1 = A*(((1-(exp(1).^(-r1*(t-t0)))).*...
%     exp(1).^(-r2*(t-t0)))./(1+((1-(exp(1).^(-r1*(t-t0)))).*...
%     exp(1).^(-r2*(t-t0)))/M)).^n.*(t>t0);
% 
% plot(t,Hog1)

% % Plot ODE solutions for mRNA:
% STL1_4state_Extended.plotODE(STL1_4state_Extended.species(5),...
%     STL1_4state_Extended.tSpan, {'linewidth',4},...
%     Title='4-state STL1 (mRNA)', TitleFontSize=24,...
%     AxisLabelSize=18, TickLabelSize=18, LegendFontSize=15,...
%     LegendLocation='east', Colors=[0.23,0.67,0.20],...
%     XLabel='Time', YLabel='Molecule Count')
% 
% STL1_4state.plotODE(STL1_4state.species(5),...
%     STL1_4state.tSpan, {'linewidth',4},...
%     Title='4-state STL1 (mRNA)', TitleFontSize=24,...
%     AxisLabelSize=18, TickLabelSize=18, LegendFontSize=15,...
%     LegendLocation='east', Colors=[0.23,0.67,0.20],...
%     XLabel='Time', YLabel='Molecule Count')
% %
% % Set 'useHybrid' to true:
% STL1_4state_Extended.useHybrid = true;
% % Define which species will be solved by ODEs:
% STL1_4state_Extended.hybridOptions.upstreamODEs = {'Hog1','D','G'};
% 
% STL1_4state_Extended.solutionScheme = 'fsp';
% STL1_4state_Extended = STL1_4state_Extended.formPropensitiesGeneral('HybridHogModel');
% [~,~,STL1_4state_Extended] = STL1_4state_Extended.solve;
% 
% % Plot the results
% STL1_4state_Extended.plotFits([], "all", [], {'linewidth',2},...
%     Title='4-state STL1', YLabel='Molecule Count',...
%     LegendLocation='northeast', LegendFontSize=12);

%% PDO
% Make FIM plots w/ w/o PDO.  
% Check if optimal expt design changes, and if so make that plot also.

% Make a copy of our model:
STL1_4state_PDO = STL1_4state_MH;

% Change to binomial PDO: 
STL1_4state_PDO_cyt = ...
    STL1_4state_PDO.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
    {'mRNA'}, {'RNA_STL1_total_TS3Full'}, {'RNA_STL1_cyto_TS3Full'},...
     'Binomial', true, [], {'Replica',1}, LegendLocation="northwest",...
     Title="4-state STL1 (Binomial PDO: Cytoplasmic mRNA)", FontSize=24,...
     XLabel="Total mRNA counts", YLabel="Cytoplasmic mRNA counts");
 
STL1_4state_PDO_nuc = ...
    STL1_4state_PDO.calibratePDO('data/filtered_data_2M_NaCl_Step.csv',...
    {'mRNA'}, {'RNA_STL1_total_TS3Full'}, {'RNA_STL1_nuc_TS3Full'},...
     'Binomial', true, [], {'Replica',1}, LegendLocation="northwest",...
     Title="4-state STL1 (Binomial PDO: Nuclear mRNA)", FontSize=24,...
     XLabel="Total mRNA counts", YLabel="Nuclear mRNA counts");

%%%%%%%%%%%%%%%%%%%%%%%%% Nucleus:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIM + PDO analyses
%   * Analyze FIM with PDO for nuclear mRNA counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Intensity
sig_log10 = 0.05*ones(1,15);
parGuess = [0, 1500, 5];

STL1_4state_PDO_intens = STL1_4state_PDO;
STL1_4state_PDO_intens = STL1_4state_PDO_intens.calibratePDO( ...
    'data/filtered_data_2M_NaCl_Step.csv', {'mRNA'},...
    {'RNA_STL1_total_TS3Full'}, {'STL1_avg_int_TS3Full'}, 'AffinePoiss',...
    true, parGuess, {'Replica',2}, LegendLocation="southeast", ...
    Title="4-state STL1 (Affine PDO: Average intensity)", FontSize=24,...
    XLabel="True mRNA counts",YLabel="Average intensities (binned)");


fimsPDOintens = STL1_4state_PDO_intens.computeFIM([],'log');
fimPDOintens = STL1_4state_PDO_intens.evaluateExperiment(fimsPDOintens,...
                                          nCellsOpt_tr,diag(sig_log10.^2));

nCellsOptPDOintens = STL1_4state_PDO_intens.optimizeCellCounts(...
                               fimsPDOintens,nTotal,'Smallest Eigenvalue');

figintens = figure;
STL1_4state_PDO_intens.plotMHResults(STL1_4state_MH_MHResults,...
                     [fimPDOintens,fimTotal,fimOpt_tr],'log',[],figintens);

%% Plot legend
axs = findall(figintens, 'Type', 'axes');
ax  = axs(1);        
hold(ax,'on');

% Helpers:
near = @(c,tol,tgt) (numel(c)==3) && all(abs(c(:)'-tgt)<=tol);
isMagenta = @(c) near(c,0.15,[1 0 1]);
isCyan    = @(c) near(c,0.15,[0 1 1]);
isBlue    = @(c) near(c,0.15,[0 0 1]);
isGreen   = @(c) near(c,0.15,[0 1 0]);

% MCMC 90% credible interval (magenta dashed):
hMHell = findobj(ax,'Type','line','LineStyle','--');
hMHell = hMHell(arrayfun(@(h) isMagenta(h.Color), hMHell));

% FIM ellipses (solid lines):
hFIM = findobj(ax,'Type','line','LineStyle','-');

% Classify FIM ellipses by color:
hFIM_cyan  = hFIM(arrayfun(@(h) isCyan(h.Color),  hFIM));
hFIM_blue  = hFIM(arrayfun(@(h) isBlue(h.Color),  hFIM));
hFIM_green = hFIM(arrayfun(@(h) isGreen(h.Color), hFIM));

% Find MCMC samples (scatter) and MLE (square marker):
hSamples = findobj(ax,'Type','scatter');
if isempty(hSamples)
    cand = findobj(ax,'Type','line','Marker','o');
    hSamples = cand(~arrayfun(@(h) strcmp(get(h,'MarkerFaceColor'),'none'), cand));
end
hMLE = findobj(ax,'Type','line','Marker','s');

% Build legend in a sensible order:
L = []; names = {};
if ~isempty(hSamples),L(end+1)=hSamples(1);names{end+1}='MH samples';end
if ~isempty(hMLE),L(end+1)=hMLE(1); names{end+1}='MLE';end
if ~isempty(hMHell),L(end+1)=hMHell(1);names{end+1}='MH 90% CI';end
if ~isempty(hFIM_cyan),L(end+1)=hFIM_cyan(1);names{end+1}='FIM PDO';end
if ~isempty(hFIM_blue),L(end+1)=hFIM_blue(1);names{end+1}='FIM total';end
if ~isempty(hFIM_green),L(end+1)=hFIM_green(1);names{end+1}='FIM optimal';end

% Fallback: if color classification failed, just take first three FIM lines
if numel(L)<5
    remainingFIM = setdiff(hFIM, [hFIM_cyan; hFIM_blue; hFIM_green]);
    for k = 1:min(3, numel(remainingFIM))
        L(end+1) = remainingFIM(k);
        names{end+1} = sprintf('FIM #%d', k);
    end
end

lgd = legend(ax, L, names, 'Location','best');
lgd.FontSize = 12;

%% Pipeline
% Brian will write a pipeline to fit 500 genes in the scSEQ data.  Jack is
% working on compiling these now to create a simplified version of existing
% data set.