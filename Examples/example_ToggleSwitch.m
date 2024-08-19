%% example_toggleSwitch
% Example script to show SSIT application to a simple two species genetic
% toggle switch.
clear all
close all
addpath(genpath('../src'));

%% Create SSIT Model
Model1 = SSIT();    % Create SSIT instance using pre-selected model
Model1.species = {'lacI';
             'lambdaCI'};   % Set species names
Model1.initialCondition = [4;0]; % Set initial condition
Model1.parameters = {'kb',10;
                'ka',80;
                'M',20;
                'g',1};   % Set parameter names and values
Model1.stoichiometry = [1,-1,0, 0;...
                   0, 0,1,-1]; % Set Stoichiometry matrix
Model1.propensityFunctions = {'kb+ka*M^3/(M^3+lambdaCI^3)';
                         'g*lacI';...
                         'kb+ka*M^3/(M^3+lacI^3)';
                         'g*lambdaCI'}; % Set propensity functions
Model1.customConstraintFuns = {'(lacI-3).^2.*(lambdaCI-3).^2'}; % Add FSP constraint
Model1.tSpan = [0:3:12];    % Set times at which to compute distributions
Model1.summarizeModel       % Print model summary to screen
Model1 = Model1.formPropensitiesGeneral('ToggleModel'); % Generate codes model

%% Solve using the FSP approach
Model1.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
Model1.fspOptions.fspTol = 1e-4;  % Set FSP error tolerance.
[Mod1FSPsoln,Model1.fspOptions.bounds] = Model1.solve;  % Solve the FSP analysis
Model1.makePlot(Mod1FSPsoln,'marginals',[2:5],false,[1,2])  % Plot marginal distributions
Model1.makePlot(Mod1FSPsoln,'joints',[2:5],false,[5])       % Plot joint distributions

%% Solve for the CME Sensitivity using FSP
Model1.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
Mod1SensSoln = Model1.solve(Mod1FSPsoln.stateSpace);  % Solve the sensitivity problem
Model1.makePlot(Mod1SensSoln,'marginals',[],false,[3,4]) % Plot marginal sensitivities

%% Define a Binomial Probabilistic Distortion Operator
Model2 = Model1;  % Make a copy of the original model
Model2.pdoOptions.type = 'Binomial';
Model2.pdoOptions.props.CaptureProbabilityS1 = 0.7;  % Distortion for S1
Model2.pdoOptions.props.CaptureProbabilityS2 = 0.7; % Distortion for S2
Model2.pdoOptions.PDO = Model2.generatePDO(Model2.pdoOptions,[0.7,0.7],Mod1FSPsoln.fsp,true);
figure(20); contourf(log10(Model2.pdoOptions.PDO.conditionalPmfs{1})); colorbar
xlabel('"true" number of mRNA'); ylabel('observed number of mRNA'); set(gca,'fontsize',15);

%% Apply PDO to FSP and Sensitivity Calculations
Model2.solutionScheme = 'FSP'; % Set solution scheme to FSP.
Model2.makePlot(Mod1FSPsoln,'marginals',[2:5],true,[1,2])  % Plot Distorted Marginals
Model2.solutionScheme = 'fspSens'; % Set solution scheme to Sensitivity
Model2.makePlot(Mod1SensSoln,'marginals',[],true,[3,4])    % Plot Distorted Sensitivities

%% Generate Simulated Data for Subsequent Fitting
Model2.solutionScheme = 'SSA';  % Set solution scheme to SSA.
Model2.ssaOptions.Nexp = 50;   % Number of independent data sets to generate.
Model2.ssaOptions.nSimsPerExpt = 200; % Number of cells to include at each time point for each data set.
Model2.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
Model2.sampleDataFromFSP(Mod1FSPsoln,'ToggleSSAData50Expts.csv')  % Generate and save all data.

%% Associate Each Dataset with an FSP Model 
Model2.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
B{1} = Model2.loadData('ToggleSSAData50Expts.csv',{'lacI','exp1_s1';'lambdaCI','exp1_s2'});
B{1}.pdoOptions.PDO = [];  % Do not use PDO.
B{2} = Model2.loadData('ToggleSSAData50Expts.csv',{'lacI','exp1_s1_Distorted';'lambdaCI','exp1_s2_Distorted'});
B{2}.pdoOptions.PDO = [];  % Do not use PDO.
B{3} = Model2.loadData('ToggleSSAData50Expts.csv',{'lacI','exp1_s1_Distorted';'lambdaCI','exp1_s2_Distorted'});

%% Sweep over Parameters to Plot Likelihood Function Landscape
fitErrorsB1 = B{1}.likelihoodSweep([2,3],linspace(.9,1.1,11),true);
title('Ideal data. Original FSP.')
% fitErrorsB2 = B2.likelihoodSweep([2,3],linspace(.5,1.5,5),true);
% title('Binomial data distortion. Original FSP.')
% fitErrorsB3 = B3.likelihoodSweep([2,3],linspace(.5,1.5,5),true);
% title('Binomial data distortion. FSP+PDO.')

%% Sweep over Parameters to Plot Likelihood Function Landscape
fitErrorsB1 = B{1}.likelihoodSweep([2,3],linspace(.5,1.5,11),true);
title('Ideal data. Original FSP.')
fitErrorsB2 = B{2}.likelihoodSweep([2,3],linspace(.5,1.5,11),true);
title('Binomial data distortion. Original FSP.')
fitErrorsB3 = B{3}.likelihoodSweep([2,3],linspace(.5,1.5,11),true);
title('Binomial data distortion. FSP+PDO.')

%% Plot Results
for m=1:3
    [fit_error,fitSolutions] = B{m}.computeLikelihood; 
    B{m}.makeFitPlot(fitSolutions); 
end

%% Find MLE for each simulated data set.
MLE = zeros(3,2,Model2.ssaOptions.Nexp); 
fMLE = inf(3,Model2.ssaOptions.Nexp);
for iExp = 1:Model2.ssaOptions.Nexp
    B{1} = B{1}.loadData('ToggleSSAData50Expts.csv',{'lacI',['exp',num2str(iExp),'_s1'];'lambdaCI',['exp',num2str(iExp),'_s2']}); % Link non-distorted data.
    B{2} = B{2}.loadData('ToggleSSAData50Expts.csv',{'lacI',['exp',num2str(iExp),'_s1_Distorted'];'lambdaCI',['exp',num2str(iExp),'_s2_Distorted']}); % Link non-distorted data.
    B{3} = B{3}.loadData('ToggleSSAData50Expts.csv',{'lacI',['exp',num2str(iExp),'_s1_Distorted'];'lambdaCI',['exp',num2str(iExp),'_s2_Distorted']}); % Link non-distorted data.
    for m=1:3
        B{m}.fittingOptions.modelVarsToFit = [2,3];
        % if iExp==1
            x0 = [B{m}.parameters{B{m}.fittingOptions.modelVarsToFit,2}]';
        % else
            % x0 = squeeze(MLE(m,:,iExp-1)); 
        % end
        [MLE(m,:,iExp),fMLE(m,iExp)] = B{m}.maximizeLikelihood(x0);
        [m;squeeze(MLE(m,:,iExp))']
    end
    % [MLE2{iExp},fMLE2(iExp)] = b2.maximizeLikelihood(x0);
    % [MLE3{iExp},fMLE3(iExp)] = b3.maximizeLikelihood(x0);
    % MLE1{iExp}
end

%% Compute FIM
cellCounts = B{1}.dataSet.nCells;  % Number of cells in each experiment.
for m=1:3
    fimResults{m} = B{m}.computeFIM; 
    % Compute the FIM for full observations and no distortion.
    [FIM{m},sFIMcov{m},fimMetrics{m}] = B{m}.evaluateExperiment(fimResults{m},cellCounts);
    fimFreePars = FIM{m}{1}(B{m}.fittingOptions.modelVarsToFit,B{m}.fittingOptions.modelVarsToFit);
    B{m}.makeMleFimPlot(squeeze(MLE(m,:,:)),fimFreePars,[1,2],0.95,100+m,[80,20])
end

%% Optimize Experiment
Model3=Model1;
Model3.tSpan = [0:1:99]; % Set times at which to compute distributions
Model3.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[Mod3SensSoln,bounds] = Model3.solve;  % Solve the sensitivity problem
cellCounts = 100*ones(size(Model3.tSpan));
nCellsTotal = sum(cellCounts);
% Compute FIM with all parameters free to vary
Model3.fittingOptions.pdoVarsToFit = 'all';
fimResults = Model3.computeFIM(Mod1SensSoln.sens); 
FIMintuit = Model3.evaluateExperiment(fimResults,cellCounts);
% Optimize cell counts for different objectives
nCellsOptDet = Model3.optimizeCellCounts(fimResults,nCellsTotal,'Determinant');
FIMoptDet= Model3.evaluateExperiment(fimResults,nCellsOptDet);
nCellsOpt1234 = Model3.optimizeCellCounts(fimResults,nCellsTotal,'[1:4]');
FIMoptDet1234= Model3.evaluateExperiment(fimResults,nCellsOpt1234);

% Plot results
pars = [[Model3.parameters{:,2}],[0.7,0.7]]';
for i=1:5
    for j = i+1:6
        subplot(5,5,(i-1)*5+j-1)
        Model3.makeMleFimPlot([],FIMintuit{1},[j,i],0.95,1,pars); hold on
        Model3.makeMleFimPlot([],FIMoptDet{1},[j,i],0.95,1,pars)
        Model3.makeMleFimPlot([],FIMoptDet1234{1},[j,i],0.95,1,pars)
        legend("off")
    end
end
legend('Intuitive','D-Optimality','$\min \det(\Sigma_{\rm model})$','interpreter','latex')
