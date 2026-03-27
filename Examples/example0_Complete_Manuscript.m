% clear

%% Instantiate SSIT
% Add SSIT source codes to Matlab search path, create an SSIT model:
addpath(genpath('../src'));
STL1_4state = SSIT('Empty');


%% 2.1.1 Define Model Species
STL1_4state.species = {'g1'; 'g2'; 'g3'; 'g4'; 'mRNA'};
 

%% 2.1.2 Define a Time-Varying Input Signal
STL1_4state.inputExpressions = ...
    {'Hog1',['A*(((1-(exp(1)^(-r1*(t-t0))))*',...
     'exp(1)^(-r2*(t-t0)))/(1+((1-(exp(1)^(-r1*(t-t0))))*',...
     'exp(1)^(-r2*(t-t0)))/M))^n*(t>t0)']};

     
%% 2.1.3 Define Reactions (Propensity Functions and Stoichiometry Matrix)
STL1_4state.propensityFunctions = {...
          'k12*g1';'(max(0,k21o*(1-k21i*Hog1)))*g2';...
          'k23*g2';'k32*g3'; 'k34*g3';'k43*g4';...
          'kr1*g1';'kr2*g2';'kr3*g3';'kr4*g4'}; 

STL1_4state.stoichiometry = [-1,1,0,0,0,0,0,0,0,0;...   % gene state 1
                              1,-1,-1,1,0,0,0,0,0,0;... % gene state 2
                              0,0,1,-1,-1,1,0,0,0,0;... % gene state 3        
                              0,0,0,0,1,-1,0,0,0,0;...  % gene state 4
                              0,0,0,0,0,0,1,1,1,1]     % mRNA
                 % Reactions: 1,2,3,4,5,6,7,8,9,10

% How to add single reaction:                 
newReaction.stoichiometry = {'mRNA',-1};  % (Reaction: 11)
newReaction.propensity = 'dr*mRNA';
newReaction.parameters = {'dr',1};
STL1_4state = STL1_4state.addReaction(newReaction);


%% 2.1.4 Define Parameters (Reaction Rates)
STL1_4state.parameters = ({'t0',3.17; 'k12',78; 'k21o',1.92e+05;...
    'k21i',3200; 'k23',0.402; 'k34',7.8; 'k32',1.62;...
    'k43',2.28; 'dr',0.294; 'kr1',4.68e-02; 'kr2',0.72;...
    'kr3', 59.4; 'kr4', 3.24; 'r1',4.14e-03; 'r2',0.426;...
    'A',9.3e+09; 'M',6.4e-04; 'n',3.1});


%% 2.1.5 Print Model Summary
STL1_4state.summarizeModel


%% 2.1.6 Save and Load a Model
save('example_1_CreateSSITModels','STL1_4state')
load('example_1_CreateSSITModels.mat')


%% 2.2 Set up for Solving Models
STL1_4state.initialCondition = [1;0;0;0;0];
STL1_4state.tSpan = linspace(0,50,101);


%% 2.2.2 Ordinary Differential Equations (ODEs)

% Set solution scheme to 'ODE':
STL1_4state.solutionScheme = 'ODE';

% Compile and store the given reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_ODE');

% Solve ODEs:
STL1_4state.Solutions = STL1_4state.solve;

% Plot ODE solutions for mRNA:
    STL1_4state.plotODE(speciesNames=STL1_4state.species(5),...
        timeVec=STL1_4state.tSpan, lineProps={'linewidth',4},...
        TitleFontSize=26, Title='4-state STL1 (mRNA)',...
        AxisLabelSize=20, TickLabelSize=20, LegendFontSize=20,...
        LegendLocation='east', Colors=[0.23,0.67,0.20],...
        XLabel='Time', YLabel='Molecule Count')

    % Plot ODE solutions for the four gene states:
    STL1_4state.plotODE(speciesNames=STL1_4state.species(1:4),...
        timeVec=STL1_4state.tSpan, lineProps={'linewidth',4},...
        TitleFontSize=26, Title='4-state STL1 (gene states)',...
        AxisLabelSize=20, TickLabelSize=20, LegendFontSize=20,...
        LegendLocation='east', XLabel='Time', YLabel='Molecule Count')


%% 2.2.3 Moment Closure

% Compile and store the given reaction propensities:
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_moments');

% Set solution scheme and solve:
STL1_4state.solutionScheme = 'moments';
[~,~,STL1_4state] = STL1_4state.solve;

% Plot moment solutions:
STL1_4state.plotMoments(solution=STL1_4state.Solutions.moments,...
    speciesNames=STL1_4state.species(5), plotType="meansanddevs",...
    lineProps={'linewidth',4}, Title='4-state STL1 (mRNA)',...
    TitleFontSize=24, YLabel='Molecule Count',...
    LegendLocation='northeast', Colors=[0.23,0.67,0.20])

STL1_4state.plotMoments(solution=STL1_4state.Solutions.moments,...
    speciesNames=STL1_4state.species(1:4), plotType="meansanddevs",...
    lineProps={'linewidth',4}, Title='4-state STL1 (gene states)',...
    TitleFontSize=24, LegendLocation='northeast', YLabel='Molecule Count')


%% 2.2.4 Stochastic Simulation Algorithm (SSA)

% Set solution scheme to SSA:
STL1_4state.solutionScheme = 'SSA';

% Set the number of simulations performed per experiment:
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
STL1_4state.Solutions = STL1_4state.solve;

% Plot SSA trajectories and means (mRNA):
STL1_4state.plotSSA(speciesIdx='all', numTraj=100,...
	speciesNames=STL1_4state.species(5), MeanOnly=false,...
    lineProps={'linewidth',4}, Title="4-state STL1 (mRNA)",...
    TitleFontSize=26, AxisLabelSize=20, TickLabelSize=20,...
    LegendFontSize=20, LegendLocation='northeast', HistTime=20,...
    XLabel='Time', YLabel='Molecule Count', Colors=[0.23,0.67,0.20]); 


%% 2.2.5 Finite State Projection (FSP)

% Ensure the solution scheme is set to FSP (default):
STL1_4state.solutionScheme = 'FSP';

% Set FSP 1-norm error tolerance:
STL1_4state.fspOptions.fspTol = 1e-4;

% Guess initial bounds on FSP StateSpace:
STL1_4state.fspOptions.bounds = [0,0,0,0,0,1,1,1,1,200];

% Approximate the steady state for the initial distribution:
STL1_4state.fspOptions.initApproxSS = true;
STL1_4state.initialTime = 0;

% Compile and store the given reaction propensities: 
STL1_4state = STL1_4state.formPropensitiesGeneral('STL1_4state_FSP');

% Solve Model:
[~,~,STL1_4state] = STL1_4state.solve;

% Plot means and standard deviations:
STL1_4state.plotFSP(speciesNames=STL1_4state.species(5),...
    plotType='meansAndDevs', lineProps={'linewidth',4},...
    Title='4-state STL1 (mRNA)', TitleFontSize=26, XLabel='Time',...
    Colors=[0.23,0.67,0.2], AxisLabelSize=20, TickLabelSize=20,...
    YLabel='Molecule Count',LegendFontSize=20,LegendLocation='northeast');

%Plot marginal distributions:
STL1_4state.plotFSP(speciesNames=STL1_4state.species(5),...
    plotType='marginals', indTimes=[1,12,24,50,101],...
    lineProps={'linewidth',3}, Colors=[0.23,0.67,0.2], XLim=[0,100])


%% 2.2.6 Escape (First Passage) and Waiting Times

% Make copy of original model:
STL1_4state_escape = STL1_4state;

% Set the initial populations:
STL1_4state_escape.initialCondition = [1;0;0;0;0];

% Turn off steady state initial condition:
STL1_4state_escape.fspOptions.initApproxSS = false;

% Set the times at which distributions will be computed:
STL1_4state_escape.tSpan = linspace(0,100,200);
STL1_4state_escape.initialTime = 0;

% Solve for time to reach mRNA=100:
STL1_4state_escape.fspOptions.escapeSinks.f = {'mRNA'};
STL1_4state_escape.fspOptions.escapeSinks.b = 100;
[~,~,STL1_4state_escape] = STL1_4state_escape.solve;

% Plot the CDF and PDF:
STL1_4state_escape.plotFSP(plotType="escapeTimes", XLim=[0,50],...
    lineProps={'linewidth',3}, LegendLocation="southeast",...
    TitleFontSize=24, Title="4-state STL1 (mRNA)", Colors=[0.23,0.67,0.2]);


%% 2.2.7 Solve FSP sensitivities

% Set solution scheme to FSP sensitivity:
STL1_4state.solutionScheme = 'fspSens';


% Solve the sensitivity problem:
[~,~,STL1_4state] = STL1_4state.solve;


% Plot the results from the sensitivity analysis
STL1_4state.plotFSP(speciesNames=STL1_4state.species(5),...
    plotType='sens', indTimes=40, lineProps={'linewidth',3},...
    Colors=[0.23,0.67,0.2], AxisLabelSize=15, TickLabelSize=12,...
    XLim=[0,100], Title="4-state STL1 (t=25)", TitleFontSize=24)


%% 2.2.8 Fisher Information Matrix (FIM) Analysis

% Set unobservable species:
STL1_4state.pdoOptions.unobservedSpecies = '1:4';

% Compute FIMs using FSP sensitivity results:
fimResults = STL1_4state.computeFIM(STL1_4state.Solutions.sens, 'log');

% Define indices of free parameters for FIM sub matrix. Here, Hog1 input 
% parameters are experimentally known (thus fixed) and all others are free:
freePars = 1:13;

% Compute the FIM sub matrix for free parameters:
fimResults_free = STL1_4state.computeFIM([],'log',[],freePars);

% Specify how many cells are to be measured at each time:
cellCounts = 1000*ones(size(STL1_4state.tSpan));

% Evaluate the provided experiment design (in "cellCounts")
% and produce an array of FIMs (one for each parameter set):
[totalFIM, STL1_4state_mleCovEstimate, fimMetrics] = ...
STL1_4state.evaluateExperiment(fimResults, cellCounts);

[freeFIM, STL1_4state_mleCovEstimate_free, fimMetrics_free] = ...
    STL1_4state.evaluateExperiment(fimResults_free, cellCounts)

% Plot the FIMs (full):
f1 = figure(11);
f2 = figure(12);
STL1_4state.plotFIMResults(totalFIM, 'log', STL1_4state.parameters,...
    PlotEllipses=true, EllipseFigure=f1,...
    EllipsePairs=[1 6; 2 3; 4 5; 6 13], FigureHandle=f2,...
    Colors=struct('EllipseColors',[0.2 0.6 0.9],...
    'CenterSquare',[0.96,0.47,0.16]));

% Plot the FIMs (free):
f3 = figure(13);
STL1_4state.plotFIMResults(freeFIM, 'log',...
    STL1_4state.parameters(1:13), PlotEllipses=true, EllipseFigure=f1,...
    EllipsePairs=[1 6; 2 3; 4 5; 6 13],FigureHandle=f3,...
    Colors=struct('EllipseColors',[0.9 0.6 0.2],...
    'CenterSquare',[0.96,0.47,0.16]));


%% 2.2.9 Experiment Design (with various FIM Optimality Criteria)
% Find the FIM-based designs for a total of 1000 cells 

% Compute the optimal number of cells from the FIM results using different 
% design criteria:  `Trace' maximizes the trace of the FIM; 
% `D-cov' minimizes the expected determinant of MLE covariance; 
% `E-opt' maximizes the smallest e.val of the FIM; and 
% `D-opt-sub[$<i_1>,<i_2>$,...]' maximizes the determinant of the FIM for  
% the specified indices.  The latter is shown for different parameter 
% combinations, where D-opt-sub[9:13]' are the mRNA-specific parameters 
% `dr' and `kr1',`kr2',`kr3', and `kr4' (degradation and transcription 
% reactions).  All other parameters are assumed to be known and fixed.
nCol = sum(cellCounts);
nTotal = nCol(1);
nCellsOpt_Dcov = STL1_4state.optimizeCellCounts(fimResults,nTotal,'D-cov');
nCellsOpt_Trace = STL1_4state.optimizeCellCounts(fimResults,nTotal,'Trace');
nCellsOpt_Doptsub = STL1_4state.optimizeCellCounts(fimResults,nTotal,...
                                                    'D-opt-sub[1:8]');
nCellsOpt_DoptsubR = STL1_4state.optimizeCellCounts(fimResults,nTotal,...
                                                    'D-opt-sub[9:13]');
nCellsOpt_DoptsubI = STL1_4state.optimizeCellCounts(fimResults,nTotal,...
                                                    'D-opt-sub[14:18]');

% Make a bar chart to compare the different designs
% Find which x positions correspond to time=30 and time=60 for off-setting:
t = STL1_4state.tSpan;
x = 1:size(t,2);   

f = figure;
bar(x,  nCellsOpt_Trace,        0.4); hold on
bar(x,  nCellsOpt_Dcov,         0.4);
bar(x,  nCellsOpt_DoptsubI,     0.4);
bar(x+0.2,  nCellsOpt_DoptsubR, 0.4);
bar(x-0.2,  nCellsOpt_Doptsub,  0.4);

set(gca,'XTick',x,'XTickLabel',t,'FontSize',16)
title('4-state STL1 (FIM Optimal Designs)','FontSize',24)
xlabel('Time (min)','FontSize',20)
ylabel('Number of cells','FontSize',20)
legend('Trace Design','D-cov Design', 'D-opt-sub[14:18] Design',...
        'D-opt-sub[9:13] Design', 'D-opt-sub[1:8] Design',...
        'Location', 'northeast')


%% 2.3.1 Data Loading and Handling

STL1_4state =STL1_4state.loadData('data/filtered_data_2M_NaCl_Step.csv',...
                                 {'mRNA','RNA_STL1_total_TS3Full'},...
                                {'Replica',1;'Condition','0.2M_NaCl_Step'});


%% 2.3.3 Maximum Likelihood Estimation
% Maximum allowable number of iterations to fit, etc.:
fitOptions = optimset('Display','iter','MaxIter',2000);

% Define which parameters to fit (in this case, all of them):
STL1_4state.fittingOptions.modelVarsToFit = [1:13];

% Search to Find the MLE:
[~,~,~,STL1_4state] = STL1_4state.maximizeLikelihood([],fitOptions);

% Make plots of the parameter fits from the MLE:
STL1_4state.plotFits(plotType="all", lineProps={'linewidth',2},...
    TitleFontSize=24, Title='4-state STL1 (MLE)', LegendFontSize=18,...
    YLabel='Molecule Count', LegendLocation='northeast', AxisLabelSize=20);

% Note: Should see an MLE of -24775.6 at the end:


%% 2.3.4 Bayesian Inference

% Make a copy of our 4-state STL1 model:
STL1_4state_MH = STL1_4state;

% Specify how many model parameters will be fit (the rest will be fixed):
fitpars = 13;

% Choose parameters to search:
STL1_4state_MH.fittingOptions.modelVarsToFit = [1:fitpars]; 

% Specify Bayesian Prior and fit:
% Specify Prior as log-normal distribution with wide uncertainty

% Prior log-mean:
mu_log10 = [0.5,2,5,3.5,-0.4,1,0.2,0.4,-0.5,-1.3,-0.1,2,0.5]; 

% Prior log-standard deviation:
sig_log10 = 2*ones(1,fitpars);  

% Prior:
STL1_4state_MH.fittingOptions.logPrior = ...
    @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Create first parameter guess:
STL1_4state_MH_pars = [STL1_4state_MH.parameters{:,2}];      

% Iterating between MLE and MH
%  Running a few rounds of MLE and MH together may improve convergence:

STL1_4state_MH.parameters(:,2) = num2cell(STL1_4state_MH_pars);

for i=1:2
    % Maximize likelihood:
    STL1_4state_MH_pars = STL1_4state_MH.maximizeLikelihood([]);    

    % Update parameters in the model:
    STL1_4state_MH.parameters([1:fitpars],2) = num2cell(STL1_4state_MH_pars);

    % Run Metropolis-Hastings    
    proposalWidthScale = 0.01;
    MHOptions.proposalDistribution  = ...
       @(x)x+proposalWidthScale*randn(size(x));

    % Set MH runtime options (number of samples, burnin, thin, etc.):
    MHOptions.numberOfSamples = 2000;
    MHOptions.burnin = 500;
    MHOptions.thin = 2;

    % Run Metropolis-Hastings: 
    [STL1_4state_MH_pars,~,STL1_4state_MH_MHResults] = ...
      STL1_4state_MH.maximizeLikelihood([],MHOptions,'MetropolisHastings');
    
    % Store MH parameters in model:
    STL1_4state_MH.parameters([1:fitpars],2) = num2cell(STL1_4state_MH_pars);
end

% Plot results:
STL1_4state_MH.plotMHResults(STL1_4state_MH_MHResults);

STL1_4state_MH.plotFits(plotType="all", lineProps={'linewidth',2},...
    Title='4-state STL1 (MH)', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=18, ProbXLim = [0 80],...
    TimePoints=[0 8 10 15 30 55], TitleFontSize=24, AxisLabelSize=20);


%% 2.3.5 Approximate Bayesian Computation

% Create a copy of our model:
STL1_4state_ABC = STL1_4state;

% Set up a prior over parameters (logPriorLoss)
logPriorLoss = @(x)sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

% Choose loss function for ABC (default: 'cdf_one_norm'):
lossFunction = 'cdf_one_norm';

% Set ABC / MCMC options
ABCoptions = struct('numberOfSamples',2000,'burnIn',500,'thin',2,...
    'proposalDistribution',@(x)x+0.03*randn(size(x)));

% Compile and store reaction propensities:
STL1_4state_ABC = STL1_4state_ABC.formPropensitiesGeneral('STL1_4state_ABC');

% Run ABC search
[~, ~, ~, STL1_4state_ABC] = STL1_4state_ABC.runABCsearch([],...
    lossFunction, logPriorLoss, ABCoptions);

STL1_4state_ABC.plotABC(STL1_4state_ABC.Solutions.ABC);


%% 2.3.6 Cross Validation

% Set Fitting Options:
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display','final','MaxIter',200); 
% Note: 'MaxIter', 200 for fast run; Set to 'MaxIter', 2000 for accuracy

% Make a copy of our 4-state STL1 model:
STL1_4state_CrossVal = STL1_4state_MH;

% Specify datafile name and species linking rules:
DataFileName = 'data/filtered_data_2M_NaCl_Step.csv';
LinkedSpecies = {'mRNA','RNA_STL1_total_TS3Full'};

% Suppose we only wish to fit the data at times before 25 minutes.  
% Set the global conditions:
ConditionsGlobal = {[],[],'TAB.time<=25'};

% Split up the replicas to be separate:
ConditionsReplicas = {'TAB.Replica==1';'TAB.Replica==2'};

% Specify constraints on rep-to-rep parameter variations. Here, we specify 
% that there is an expected 0.1 log10 deviation expected in some parameters 
% and smaller in others.  No deviation at all is indicated by 0.
Log10Constraints = ...
    [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.02,0.02,0.02,0.02,0.02]; 

% Create full model:
CrossValidationModel = SSITMultiModel.createCrossValMultiModel(...
    STL1_4state_CrossVal, DataFileName, LinkedSpecies, ConditionsGlobal,...
    ConditionsReplicas, Log10Constraints);
CrossValidationModel = CrossValidationModel.initializeStateSpaces;

% Run the model fitting routines:
crossValPars = CrossValidationModel.parameters;
crossValPars = CrossValidationModel.maximizeLikelihood(...
    crossValPars, fitOptions, fitAlgorithm);
CrossValidationModel = CrossValidationModel.updateModels(crossValPars);
CrossValidationModel.parameters = crossValPars;

% Make a figure to explore how much the parameters changed between replicas:
fignum = 11; useRelative = true;
CrossValidationModel.compareParameters(fignum,useRelative);


%% 2.4 Prelim (Define and solve scRNA-seq template model):
% Define Base Model Combination
scRNAseq = SSIT;
scRNAseq.species = {'onGene';'rna'};
scRNAseq.initialCondition = [0;0];
scRNAseq.propensityFunctions = ...
    {'(kon_0+kon_1*Iupstream)*(2-onGene)';...
    'koff_0/(1+akoff*Iupstream)*onGene';...
    'kr_0*(2-onGene)+kr_1*onGene';'gr*rna'};
scRNAseq.stoichiometry = [1,-1,0,0;0,0,1,-1];
scRNAseq.inputExpressions = {'Iupstream',...
                                'exp(-r1*t*(t>=0))*(1-exp(-r2*t*(t>=0)))'};
scRNAseq.parameters = ({'r1',0.01; 'r2',0.1; 'kon_0',0.01;...
                              'kon_1',0.01; 'koff_0',20; 'akoff',0.2;...
                              'kr_0',1; 'kr_1',100; 'gr',1});

scRNAseq.fspOptions.initApproxSS = true;
scRNAseq.fittingOptions.modelVarsToFit = 1:9;
scRNAseq.fittingOptions.logPrior = @(x)-sum(log10(x).^2/2);

% We generate functions for model propensities
scRNAseq = scRNAseq.formPropensitiesGeneral('scRNAseq_Template');
[~,~,scRNAseq] = scRNAseq.solve;


%% 2.4.1 Model Reduction
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
STL1_ModRed.plotFSP(solution=redSoln, plotType='meansAndDevs',...
    speciesNames=STL1_ModRed.species(5), lineProps={'linewidth',4},...
    Title='4-state STL1 (FSP Reduced)', TitleFontSize=24,...
    XLabel='Time', Colors=[0.23,0.67,0.2], YLabel='Molecule Count',...
    LegendFontSize=15, LegendLocation='northeast',YLim=[0,40]);

% Plot the full and reduced FSP solutions:
STL1_4state.plotFSP(solution=fullSoln, plotType='meansAndDevs',...
    speciesNames=STL1_ModRed.species(5), lineProps={'linewidth',4},...
    Title='4-state STL1 (FSP Full)', TitleFontSize=24,...
    XLabel='Time', Colors=[0.23,0.67,0.2], YLabel='Molecule Count',...
    LegendFontSize=15, LegendLocation='northeast',YLim=[0,40]);


%% 2.4.2 Hybrid Models
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
STL1_4state_Extended.plotFits(plotType="all", lineProps={'linewidth',2},...
    Title='4-state STL1', YLabel='Molecule Count',...
    LegendLocation='northeast', LegendFontSize=12);


%% 2.4.3 Data Distortion Handling (PDOs)
% Handle data distortions with probability distortion operators (PDOs)

% Set PDO to Binomial for RNA seq data (assuming 95% drop-out rate)
scRNAseq.pdoOptions.type = 'Binomial';
scRNAseq.pdoOptions.unobservedSpecies = 'onGene';
scRNAseq.pdoOptions.props.CaptureProbabilityS1 = 0;    % Gene State is not measured
scRNAseq.pdoOptions.props.CaptureProbabilityS2 = 0.05; % 95% drop out from RNA
[~,scRNAseq] = scRNAseq.generatePDO(showPlot=true,...
    Title='RNA Seq (Binomial PDO: 95% Drop Out)');

% Make a copy of our 4-state STL1 model:
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


%% 2.4.4 Multi Model

% Generate library of individual gene models
% Specify datafile name:
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';

% Get name of each gene for species linking:
TAB = readtable(DataFileName);
geneNames = fields(TAB);
geneNames = geneNames(2:end-4); % the genes are in columns 2 -> N-4

if ~exist('seqModels','dir'); mkdir('seqModels'); end

% Link species 'rna' to RNA count for each gene:
for iGene = 1:length(geneNames)
    linkedSpecies = {'rna',geneNames{iGene}};
    Model = Model_Template.loadData(DataFileName,linkedSpecies);
    modelName = ['Model_',geneNames{iGene}];
    assignin('base',modelName,Model);
    save(['seqModels/',modelName],modelName);
end

% Fit a multi-model for the four genes (constrain parameters of the 
% upstream input signal):

% Select which models to include in multimodel.
Models = {Model_DUSP1,Model_RUNX1,Model_BIRC3,Model_TSC22D3};
modelNames = {'Model_DUSP1','Model_RUNX1','Model_BIRC3','Model_TSC22D3'};

% Define how parameters are assigned to sub-models (all genes are
% assumed to use the same upstream signal, but have different gene bursting
% parameters):
ParInds = {[1,2,3:9],[1,2,7*1+(3:9)],[1,2,7*2+(3:9)],[1,2,7*3+(3:9)]};

% Define constraint on model parameters (should be of similar
% magnitudes for each gene unless otherwise demanded by the differences in
% the data):
Constraint = @(x) -var(log10([x(3:9);x(7*1+(3:9));x(7*2+(3:9));x(7*3+(3:9))]));

% Create and initialize multimodel
combinedModel = SSITMultiModel(Models, ParInds, Constraint);
combinedModel = combinedModel.initializeStateSpaces();

% Fit the multimodel:
% Because there are a lot of parameters, this could take a few rounds to
% get a good MLE. You should be able to get a best posterior of about XXX.
fitOptions = optimset('Display','iter','MaxIter',1000);
fitOptions.suppressExpansion = true;
for i = 1:10
    [~,~,~,combinedModel] = combinedModel.maximizeLikelihood([],fitOptions);
    save('seqModels/CombinedModel4Genes','combinedModel');
end


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