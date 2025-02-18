%% example_DUSP1_Regulation
% Example script to show SSIT application to a simple bursting gene
% expression model that ios being fit to smFISH data for DUSP1 upon
% stimulation with Dexamethsone.
clear; clc; close all
addpath(genpath('../src'));

%% STEP1 == Define SSIT Model
% Here we set up a simple model where there is an upstream transcription
% factor (GR) that activates a gene.  Once active, the gene can transcribe
% nuclear RNA, which can later decay or leave the nucleus.
Model1 = SSIT;  % Create blank SSIT model.
Model1.species = {'offGene';'onGene';'rna'}; % Set species names.
Model1.initialCondition = [2;0;0];           % Set Initial condition

% Define propensity functions and input signals:
Model1.propensityFunctions = {'kon*IGR*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
Model1.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'}; 
Model1.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; % Define stoichiometry
Model1.parameters = ({'koff',0.014;'kon',0.002;'kr',1;'gr',0.004;...
                      'a1',20;'r1',0.04;'r2',0.1});  % Specify parameter guesses
Model1.fspOptions.initApproxSS = true;  % Set Initial Distribution to Steady State.
Model1.summarizeModel                   % Print visual summary of model
Model1 = Model1.formPropensitiesGeneral('ToyDUSP1Model'); % Generate model codes

%% STEP2 == Solve CME using FSP
% Next, we can solve the model using the FSP.  In this example, we show how
% to run the code twice.  First call finds the FSP projection needed to
% solve the problem, and the second call solves using that projection.
Model1.solutionScheme = 'FSP';   % Select FSP solution Scheme
Model1.fspOptions.fspTol = 1e-4; % Set FSP 1-norm error tolerance.
Model1.fspOptions.bounds(4:6) = [2,2,400];  % Guess initial bounds on FSP StateSpace
Model1.tSpan = linspace(0,180,301);
[Mod1FSPsoln,Model1.fspOptions.bounds] = Model1.solve; % Solve Model
Model1.makePlot(Mod1FSPsoln,'marginals',[1:100:301],false,[1,2,3],{'linewidth',2})  % Plot marginal distributions
Model1.makePlot(Mod1FSPsoln,'margmovie',[],false,[101],{'linewidth',2},'movie.mp4',[1,1,0.015],[2,3])  % Plot marginal distributions

%% STEP3 == Solve Sensitivity using FSP
Model1.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
Mod1SensSoln = Model1.solve(Mod1FSPsoln.stateSpace);  % Solve the sensitivity problem
Model1.makePlot(Mod1SensSoln,'marginals',[1:100:301],false,[3+(1:12)],{'linewidth',2}) % Plot marginal sensitivities
Model1.makePlot(Mod1SensSoln,'margmovie',[],false,[101],{'linewidth',2},'sensMovie.mp4',[-1,-1,-0.015;1,1,0.015],[3],[3])  % Plot marginal distributions

%% STEP4 == Define Binomial Probabilistic Distortion Operator
Model2 = Model1;  % Make a copy of the original model
Model2.pdoOptions.type = 'Binomial';
Model2.pdoOptions.props.CaptureProbabilityS1 = 0;  % Distortion for OFF species (unobserved)
Model2.pdoOptions.props.CaptureProbabilityS2 = 0;  % Distortion for ON species (unobserved)
Model2.pdoOptions.props.CaptureProbabilityS3 = 0.7;% Distortion for RNA species
Model2.pdoOptions.PDO = Model2.generatePDO(Model2.pdoOptions,[],Mod1SensSoln.sens.data,true);
figure(20); contourf(log10(Model2.pdoOptions.PDO.conditionalPmfs{3}),30); colorbar
xlabel('"true" number of mRNA'); ylabel('observed number of mRNA'); set(gca,'fontsize',15);

%% STEP5 == Apply PDO to FSP and Sensitivity Calculations
Model2.solutionScheme = 'FSP'; % Set solution scheme to FSP.
Model2.makePlot(Mod1FSPsoln,'marginals',[1:100:301],true,[1,2,3],{'linewidth',2})  % Plot Distorted Marginals
Model2.solutionScheme = 'fspSens'; % Set solution scheme to Sensitivity
Model2.makePlot(Mod1SensSoln,'marginals',[1:100:301],true,[3+(1:12)],{'linewidth',2})    % Plot Distorted Sensitivities

%% STEP6 == Generate Simulated Data for Fitting
Model2.tSpan = [0,60,120,180]; % Set the times at which to generate data.
Model2.solutionScheme = 'FSP'; % Set solution scheme to FSP.
Mod2FSPsoln = Model2.solve; % Solve Model
Model2.ssaOptions.Nexp = 50;   % Number of independent data sets to generate.
Model2.ssaOptions.nSimsPerExpt = 200; % Number of cells to include at each time point for each data set.
Model2.ssaOptions.applyPDO = true;    % Include the distortion in the data.
dataTable = Model2.sampleDataFromFSP(Mod2FSPsoln,'DUSP1SSAData50Expts.csv'); % Generate and save data.

% Plot data as histograms
tPlot = [0,60,120,180];
for i = 1:4
    subplot(2,2,i)  % Switch to current subplot.
    histogram(dataTable.exp1_s3(dataTable.time==tPlot(i)),30,"DisplayStyle","stairs"); hold on;
    histogram(dataTable.exp1_s3_Distorted(dataTable.time==tPlot(i)),30,"DisplayStyle","stairs");
end

%% STEP7 == Associate Datasets with FSP Models
Model2.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
B{1} = Model2.loadData('DUSP1SSAData50Expts.csv',{'rna','exp1_s3'});
B{1}.pdoOptions.PDO = [];  % Do not use PDO.
B{2} = Model2.loadData('DUSP1SSAData50Expts.csv',{'rna','exp1_s3_Distorted'});
B{2}.pdoOptions.PDO = [];  % Do not use PDO.
B{3} = Model2.loadData('DUSP1SSAData50Expts.csv',{'rna','exp1_s3_Distorted'});

%% STEP8 == Sweep overParameters and Plot Likelihood Functions
fitErrorsB1 = B{1}.likelihoodSweep([2,3],linspace(.5,1.5,11),true);
title('Ideal data. Original FSP.')
fitErrorsB2 = B{2}.likelihoodSweep([2,3],linspace(.5,1.5,11),true);
title('Binomial data distortion. Original FSP.')
fitErrorsB3 = B{3}.likelihoodSweep([2,3],linspace(.5,1.5,11),true);
title('Binomial data distortion. FSP+PDO.')

%% STEP9 == Find MLEs for Simulated Datasets
for iExp = 1:Model2.ssaOptions.Nexp
    for m = 1:3
        switch m % Link appropriate data sets
            case 1; b = B{1}.loadData('DUSP1SSAData50Expts.csv',{'rna',['exp',num2str(iExp),'_s3']}); 
            case 2; b = B{2}.loadData('DUSP1SSAData50Expts.csv',{'rna',['exp',num2str(iExp),'_s3_Distorted']}); 
            case 3; b = B{3}.loadData('DUSP1SSAData50Expts.csv',{'rna',['exp',num2str(iExp),'_s3_Distorted']});
        end
        b.fittingOptions.modelVarsToFit = [2,3];   % Only fit two of the parameters
        x0 = [b.parameters{b.fittingOptions.modelVarsToFit,2}]';  % Initial parameter guess
        [MLE(m,:,iExp),fMLE(m,iExp)] = b.maximizeLikelihood(x0);  % Run fitting code
    end
end

%% STEP10 == Compute FIM and Compare to MLE Spread
for m=1:3
    fimResults{m} = B{m}.computeFIM;  % Compute FIM for individual times
    % Compute the FIM for full observations
    FIM{m} = B{m}.evaluateExperiment(fimResults{m},B{m}.dataSet.nCells);
    iVars = b.fittingOptions.modelVarsToFit; % Indices of free variables.
    fimFreePars = FIM{m}{1}(iVars,iVars);    % Select FIM for parameters of interest
    b.makeMleFimPlot(squeeze(MLE(m,:,:)),fimFreePars,[1,2],.95,100+m,[0.002,1]) % Plot results
end

%% STEP11 == Add real data to the model
mReal = Model1;  % Make copy of previous model
mReal = mReal.loadData('../ExampleData/Dusp1Data.csv',...
    {'rna','RNA_DUSP1_nuc'},{'Dex_Conc','100'}); % Load data into model

%% STEP12a == Compute FIM, Run Metropolis Hastings
% Specify Prior as log-normal distribution with wide uncertainty
mu_log10 = [-1,-2,0,-2,1,-1,-1];  % Prior log-mean
sig_log10 = 2*ones(1,7);          % Prior log-standard deviation
mReal.fittingOptions.logPrior = @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

mReal.fittingOptions.modelVarsToFit = [1:7]; % Choose parameters to search
DUSP1pars = [mReal.parameters{:,2}];         % Create first parameter guess


fimResults = mReal.computeFIM([],'log'); % Compute individual FIMs
fimTotal = mReal.evaluateExperiment(fimResults,mReal.dataSet.nCells,...
    diag(sig_log10.^2)); % Compute total FIM including effect of prior.
mReal.fittingOptions.modelVarsToFit = [1:4]; % Choose parameters to search
FIMfree = fimTotal{1}([1:4],[1:4]); % Select FIM for free parameters.
COVfree = (1/2*(FIMfree+FIMfree'))^(-1);  % Estimate Covariance using CRLB.

% Define Metropolis Hasting Settings.
mReal.fittingOptions.logPrior = @(x)-sum((log10(x)-mu_log10([1:4])).^2./(2*sig_log10([1:4]).^2));
MHFitOptions = struct('proposalDistribution',@(x)mvnrnd(x,COVfree),...
    'numberOfSamples',5000);
[DUSP1pars,~,MHResultsDusp1] = mReal.maximizeLikelihood(...
    [], MHFitOptions, 'MetropolisHastings'); % Run Metropolis Hastings
mReal.parameters([1:4],2) = num2cell(DUSP1pars);

mReal.plotMHResults(MHResultsDusp1,FIMfree,'log',[])

%% STEP12b == Specify Bayesian Prior and fit.
% Specify Prior as log-normal distribution with wide uncertainty
mu_log10 = [-1,-2,0,-2,1,-1,-1];  % Prior log-mean
sig_log10 = 2*ones(1,7);          % Prior log-standard deviation
mReal.fittingOptions.logPrior = @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));

mReal.fittingOptions.modelVarsToFit = [1:7]; % Choose parameters to search
DUSP1pars = [mReal.parameters{:,2}];         % Create first parameter guess
DUSP1pars = mReal.maximizeLikelihood(DUSP1pars); % Fit to maximize likelihood
mReal.parameters(:,2) = num2cell(DUSP1pars); % Update new parameters.
mReal.makeFitPlot  % Plot fitting results
% You may need to re-run this multiple times until converged.
% I got a MLE of -52,454.1 after a few runs. 

%% STEP13 == Compute FIM, Run Metropolis Hastings
fimResults = mReal.computeFIM([],'log'); % Compute individual FIMs
fimTotal = mReal.evaluateExperiment(fimResults,mReal.dataSet.nCells,...
    diag(sig_log10.^2)); % Compute total FIM including effect of prior.
mReal.fittingOptions.modelVarsToFit = [1:4]; % Choose parameters to search
FIMfree = fimTotal{1}([1:4],[1:4]); % Select FIM for free parameters.
COVfree = (1/2*(FIMfree+FIMfree'))^(-1);  % Estimate Covariance using CRLB.

% Define Metropolis Hasting Settings.
mReal.fittingOptions.logPrior = @(x)-sum((log10(x)-mu_log10([1:4])).^2./(2*sig_log10([1:4]).^2));
MHFitOptions = struct('proposalDistribution',@(x)mvnrnd(x,COVfree),...
    'numberOfSamples',5000);
[DUSP1pars,~,MHResultsDusp1] = mReal.maximizeLikelihood(...
    [], MHFitOptions, 'MetropolisHastings'); % Run Metropolis Hastings
mReal.parameters([1:4],2) = num2cell(DUSP1pars);

mReal.plotMHResults(MHResultsDusp1,FIMfree,'log',[])

%% STEP14 == Optimize Experiment 
nCellsOriginal = mReal.dataSet.nCells; % Original Expt Design
nCellsTotal = sum(nCellsOriginal); % Compute total number of cells
% Optimize cell counts to minimize uncertainty in parameters 1-4.
nCellsOpt = mReal.optimizeCellCounts(fimResults,nCellsTotal,'Determinant',...
    [],[],[],[],diag(sig_log10.^2));

% Plot experiment designs
figure()
bar([1:12],nCellsOriginal,0.4); hold on;
bar([1:12]+0.5,nCellsOpt,0.4); hold on;
set(gca,'fontsize',15,'xtick',[1:12]+.25,'xticklabel',mReal.tSpan)

%% STEP15 == Fit reduced data set
mRealReduced = mReal;   % Copy previous model
mRealReduced.fittingOptions.timesToFit = nCellsOpt>0; % Ignore discarded timepoints.
mu_log10 = [-1,-2,0,-2,1,-1,-1];  % Prior log-mean
sig_log10 = 2*ones(1,7);          % Prior log-standard deviation
mRealReduced.fittingOptions.logPrior = @(x)-sum((log10(x)-mu_log10).^2./(2*sig_log10.^2));
parsRed = [mRealReduced.parameters{:,2}];  % Create first parameter guess
parsRed = mRealReduced.maximizeLikelihood(parsRed); % Fit to maximize likelihood
mRealReduced.parameters(:,2) = num2cell(parsRed);   % Update new parameters.

mRealReduced.fittingOptions.timesToFit = ones(size(nCellsOpt),'logical');
%%
%mRealReduced.makeFitPlot()  % Plot fitting results

%% STEP16 == Calibrate PDO from Eric's DUSP1 Intensity Data
ModelPDOIntensEric = mReal;
ModelPDOIntensEric = ModelPDOIntensEric.calibratePDO('../ExampleData/Dusp1Data.csv',...
    {'rna'},{'RNA_DUSP1_nuc'},{'Nuc_DUSP1_avg_int_tot'},'AffinePoiss',true,[1,230,0.5]);
