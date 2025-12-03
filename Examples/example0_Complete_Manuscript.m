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

% Solve ODEs:
[~,~,STL1_4state] = STL1_4state.solve;

% Plot ODE solutions for mRNA:
STL1_4state.plotODE(STL1_4state.species(5),...
    STL1_4state.tSpan, {'linewidth',4},...
    Title='4-state STL1 (mRNA)', TitleFontSize=24,...
    AxisLabelSize=18, TickLabelSize=18, LegendFontSize=15,...
    LegendLocation='east', Colors=[0.23,0.67,0.20],...
    XLabel='Time', YLabel='Molecule Count',...
    XLim=[0,50], YLim=[0,25])

% Plot ODE solutions for the 4 gene states:
STL1_4state.plotODE(STL1_4state.species(1:4),...
    STL1_4state.tSpan, {'linewidth',4},...
    Title='4-state STL1 (gene states)', TitleFontSize=24,...
    AxisLabelSize=18, TickLabelSize=18,...
    LegendFontSize=15, LegendLocation='east',...
    XLabel='Time', YLabel='Molecule Count',...
    XLim=[0,50], YLim=[0,1])

%% 2.10
% Set solution scheme and solve.
STL1_4state.solutionScheme = 'moments';
[~,~,STL1_4state] = STL1_4state.solve;

% TODO - Alex, please make a new plotMoments function with the following,
% but cleaned up like the other plots you have generated:

% Number of species
nSp = numel(STL1_4state.species);

% Extract mean and second moment for STL1
STL1_means = STL1_4state.Solutions.moments(nSp, :);
STL1_second = STL1_4state.Solutions.moments(end, :);

% Compute variance of STL1
STL1_var =STL1_second - STL1_means.^2;

figure
errorbar(STL1_4state.tSpan,STL1_means,sqrt(STL1_var))

% TODO - Make plots of means +/- STD vs time using results.

%% 2.11
% Set solution scheme to SSA:
STL1_4state.solutionScheme = 'SSA';

% Set the number of simulations performed per experiment
STL1_4state.ssaOptions.Nsims = 100;

% A negative initial time can be used to allow model to equilibrate (burn-in)
% to an initial steady state prior to starting the subsequent simulation.
STL1_4state.tSpan = [-100,STL1_4state.tSpan];

% Set the initial time:
STL1_4state.initialTime = STL1_4state.tSpan(1);

% Run iterations in parallel with multiple cores, or execute serially:
STL1_4state.ssaOptions.useParallel = true;

% Run SSA:
[~,~,STL1_4state] = STL1_4state.solve;

% Plot SSA trajectories and means (mRNA):
STL1_4state.plotSSA('all',100,STL1_4state.species(5),{'linewidth',4}, HistTime=20,...
Title="4-state STL1 (mRNA)", MeanOnly=false,...
Colors=[0.23,0.67,0.20], TitleFontSize=24, LegendLocation='northeast');

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

% Solve Model:
[~,~,STL1_4state] = STL1_4state.solve;

% Plot means and standard deviations:
STL1_4state.plotFSP(STL1_4state.Solutions,...
    STL1_4state.species(5), 'meansAndDevs', [], [], {'linewidth',4},...
    Title='4-state STL1 (mRNA)', TitleFontSize=24, LegendFontSize=15,...
    Colors=[0.23,0.67,0.2], LegendLocation='northeast');

% Plot marginal distributions:
STL1_4state.plotFSP(STL1_4state.Solutions,...
    STL1_4state.species(5), 'marginals', [1,2,3,5,8,13],...
    [], {'linewidth',3}, Colors=[0.23,0.67,0.2], XLim=[0,100])
%% FSP Escape Times
% Make Copy of Original Model
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

% Solve the sensitivity problem:
[~,~,STL1_4state] = STL1_4state.solve;

% Choose time at which to make sensitivity plots
index_PlotTime = 6;  

% Plot the results from the sensitivity analysis
STL1_4state.plotFSP(STL1_4state.Solutions,...
    STL1_4state.species(5), 'sens', index_PlotTime, [], {'linewidth',3}, ...
    Colors=[0.23,0.67,0.2], AxisLabelSize=14, TickLabelSize=10, ...
    XLim=[0,10], Title="4-state STL1 (mRNA)", TitleFontSize=18)

%% FIM Analysis

% Set unobservable species.

% Add lognormal prior (so FIM is invertible).

% Compute FIMs using FSP sensitivity results.
fimResults = STL1_4state.computeFIM;

% Specify how many cells are to be measured at each time:
STL1_4state_cellCounts = 1000*ones(size(STL1_4state.tSpan));

% Evaluate the provided experiment design (in "cellCounts")
% and produce an array of FIMs (one for each parameter set):
[STL1_4stateTotalFIM,...
    STL1_4state_mleCovEstimate,STL1_4stateMetrics] = ...
    STL1_4state.evaluateExperiment(fimResults,...
    STL1_4state_cellCounts);

% Plot the FIMs:
STL1_4state.plotFIMResults(STL1_4stateTotalFIM,...
    STL1_4state.parameters, [STL1_4state.parameters{:,2}],...
    PlotEllipses=true, EllipsePairs=[1 9; 3 7; 8 10; 6 9]);

% TODO - Alex, please note that the order of the parameters has changed
% from the figutre you had before -- this reflects the order in which the
% parameters were originally defined in Box 2.5.


%% Experiment Design

% Choose three design criteria (D, E, Ds) -- find the three designs for
% totla of 1000 cells -- make plots of expt designs in bar charts (like
% that already in paper).

%% Load and Plot Data
STL1_4state.solutionScheme = 'fsp';

% TODO - Make sure to add a comment saying that the user needs to change the path
% to match where the data is.
STL1_4state = ...
    STL1_4state.loadData('data/filtered_data_2M_NaCl_Step.csv',...
    {'mRNA','RNA_STL1_total_TS3Full'},...
    {'Replica',2;'Condition','0.2M_NaCl_Step'});
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

% TODO - Add comment recording the lowest found MLE.

%%
%% Specify Bayesian Prior and fit
% Specify Prior as log-normal distribution with wide uncertainty
% Prior log-mean:
mu_log10 = [0.8,3,-0.1,2,2.75,0.6,3,2.5,0,3.5,1.5,-0.15,0.5,1.5,-1];

% Prior log-standard deviation:
sig_log10 = 2*ones(1,15);

% Prior:
STL1_4state.fittingOptions.logPrior = ...
@(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));
STL1_4state.fittingOptions.logPriorCovariance = diag(sig_log10);

% Fit to maximize likelihood:
[~,~,~,STL1_4state] = STL1_4state.maximizeLikelihood;

% Plot fitting results:
STL1_4state.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);

% Iterating between MLE and MH
% Running a few rounds of MLE and MH together may improve convergence.

for i=1:3
    % Maximize likelihood:
    [~,~,~,STL1_4state] = STL1_4state.maximizeLikelihood;
    
    % Set MH runtime options (number of samples, burnin, thin, etc.):
    MHOptions = struct('numberOfSamples',2000,...
    'burnin',200,'thin',2,'useFIMforMetHast',true);
   
    % Run Metropolis-Hastings:    
    [~,~,~,STL1_4state] = ...
        STL1_4state.maximizeLikelihood([], MHOptions,...
        'MetropolisHastings');
end
STL1_4state.plotMHResults;
STL1_4state.plotFits([], "all", [], {'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12, ProbXLim = [0 80]);

% TODO - Need a wrapper for the MH Plotting functions to allow easier
% choice of which scatter plots to show.

%% ABC
STL1_4state_ABC = STL1_4state;
% Set up a prior over parameters (logPriorLoss)
logPriorLoss = @(x)sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));
% Choose loss function for ABC (default: 'cdf_one_norm'):
lossFunction = 'cdf_one_norm';

% Set ABC / MCMC options
ABCoptions = struct('numberOfSamples',2000,'burnIn',0,'thin',1,...
    'proposalDistribution',@(x)x+0.05*randn(size(x)));

% Run ABC search
[~, ~, ~, STL1_4state_ABC] = ...
    STL1_4state_ABC.runABCsearch([], lossFunction, logPriorLoss,...
    ABCoptions);

STL1_4state_ABC.plotMHResults(STL1_4state_ABC.Solutions.ABC);


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

% Create full model:
CrossValidationModel = SSITMultiModel.createCrossValMultiModel(...
    STL1_4state, DataFileName, LinkedSpecies, ConditionsGlobal,...
    ConditionsReplicas, Log10Constraints);

% Run the model fitting routines:
[~,~,~,CrossValidationModel] = CrossValidationModel.maximizeLikelihood([],...
    fitOptions);

% Make a figure to explore how much the parameters changed between replicas:
fignum = 12; useRelative = true;
CrossValidationModel.compareParameters(fignum,useRelative);


%% Multi-Model
% This looks very similar to the previous cross-validation model. I recommend
% combining it with the above.

%% Model Reduction
% None of the current model reductions are meant for use in time varying
% problems, so I doubt that they would work for the Hog Model.  I recommend
% demonstrating this on a different model.

% Make a copy of the STL1 model to set up for model reduction:
STL1_ModRed = STL1_4state;
STL1_ModRed.fspOptions.initApproxSS = true;
STL1_ModRed.modelReductionOptions.useModReduction = true;
STL1_ModRed.fspOptions.fspTol = inf;
STL1_ModRed.modelReductionOptions.reductionType = 'Logarithmic State Lumping';
STL1_ModRed.modelReductionOptions.reductionOrder = 40;
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
% Change to binomial PDO. Plot PDO, make FSP solutions plots with and
% without PDO. Make FIM plots w/ w/o PDO.  Check if optimal expt design
% changes, and if so make that plot also.

%% Pipeline
% Brian will write a pipeline to fit 500 genes in the scSEQ data.  Jack is
% working on compiling these now to create a simplified version of existing
% data set.