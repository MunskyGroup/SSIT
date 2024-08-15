%% example_toggleSwitch
% Example script to show SSIT application to a simple two species genetic
% toggle switch.
clear all
close all
addpath(genpath('../src'));

%% Create SSIT Model
A = SSIT();    % Create SSIT instance using pre-selected model
A.species = {'lacI';
             'lambdaCI'};   % Set species names.
A.parameters = {'kb',10;
                'ka',80;
                'M',20;
                'g',1}; % Set parameter names and values
A.stoichiometry = [1,-1,0, 0;...
                   0, 0,1,-1]; % Set Stoichiometry matrix
A.propensityFunctions = {'kb+ka*M^3/(M^3+lambdaCI^3)';
                         'g*lacI';...
                         'kb+ka*M^3/(M^3+lacI^3)';
                         'g*lambdaCI'}; % Set propensity functions
A.initialCondition = [4;0]; % Set initial condition
A.customConstraintFuns = {'(lacI-3).^2.*(lambdaCI-3).^2'}; % Set FSP constraint.
A.tSpan = [0:3:12];          % Set times at which to compute distributions
A = A.formPropensitiesGeneral('ToggleModel');

%% Solve using the FSP approach
A.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
A.fspOptions.fspTol = 1e-4;  % Set FSP error tolerance.
[FSPsoln,A.fspOptions.bounds] = A.solve;  % Solve the FSP analysis
A.makePlot(FSPsoln,'marginals',[2:5],false,[1,2])  % Plot marginal distributions
A.makePlot(FSPsoln,'joints',[2:5],false,[5])       % Plot joint distributions

%% Solve Sensitivity using FSP
A.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln] = A.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem
A.makePlot(sensSoln,'marginals',[],false,[3,4]) % Plot marginal sensitivities

%% Define a Binomial PDO
A.pdoOptions.type = 'Binomial';
A.pdoOptions.props.CaptureProbabilityS1 = 0.7;  % Distortion for S1
A.pdoOptions.props.CaptureProbabilityS2 = 0.7; % Distortion for S2
A.pdoOptions.PDO = A.generatePDO(A.pdoOptions,[0.7,0.7],FSPsoln.fsp,true);
figure(20); contourf(log10(A.pdoOptions.PDO.conditionalPmfs{1})); colorbar
xlabel('"true" number of mRNA'); ylabel('observed number of mRNA'); set(gca,'fontsize',15);

%% Apply PDO to FSP and Sensitivity Calculations
A.solutionScheme = 'FSP'; % Set solution scheme to FSP.
A.makePlot(FSPsoln,'marginals',[2:5],true,[1,2])  % Plot Distorted Marginals
A.solutionScheme = 'fspSens'; % Set solution scheme to Sensitivity
A.makePlot(sensSoln,'marginals',[],true,[3,4])    % Plot Distorted Sensitivities

%% Simulate Data for Subsequent Fitting
A.solutionScheme = 'FSP';  % Set solution scheme to SSA.
A.ssaOptions.Nexp = 150;   % Number of independent data sets to generate.
A.ssaOptions.nSimsPerExpt = 200; % Number of cells to include at each time point for each data set.
A.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
A.sampleDataFromFSP(FSPsoln,'ToggleSSAData50Expts.csv')  % Generate and save all data.

%% Associate Datasets with FSP Models
A.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
B{1} = A.loadData('ToggleSSAData50Expts.csv',{'lacI','exp1_s1';'lambdaCI','exp1_s2'});
B{1}.pdoOptions.PDO = [];  % Do not use PDO.
B{2} = A.loadData('ToggleSSAData50Expts.csv',{'lacI','exp1_s1_Distorted';'lambdaCI','exp1_s2_Distorted'});
B{2}.pdoOptions.PDO = [];  % Do not use PDO.
B{3} = A.loadData('ToggleSSAData50Expts.csv',{'lacI','exp1_s1_Distorted';'lambdaCI','exp1_s2_Distorted'});

%% Plot Results
for m=1:3; [fit_error,fitSolutuions] = B{m}.computeLikelihood; B{m}.makeFitPlot(fitSolutuions); end

%% Find MLE for each simulated data set.
MLE = zeros(3,2,A.ssaOptions.Nexp); fMLE = inf(3,2,A.ssaOptions.Nexp);
for iExp = 1:A.ssaOptions.Nexp
    B{1} = B{1}.loadData('ToggleSSAData50Expts.csv',{'lacI',['exp',num2str(iExp),'_s1'];'lambdaCI',['exp',num2str(iExp),'_s2']}); % Link non-distorted data.
    B{2} = B{2}.loadData('ToggleSSAData50Expts.csv',{'lacI',['exp',num2str(iExp),'_s1_Distorted'];'lambdaCI',['exp',num2str(iExp),'_s2_Distorted']}); % Link non-distorted data.
    B{3} = B{3}.loadData('ToggleSSAData50Expts.csv',{'lacI',['exp',num2str(iExp),'_s1_Distorted'];'lambdaCI',['exp',num2str(iExp),'_s2_Distorted']}); % Link non-distorted data.
    for m=1:3
        B{m}.fittingOptions.modelVarsToFit = [2,3];
        if iExp==1; x0 = [B{m}.parameters{B{m}.fittingOptions.modelVarsToFit,2}]';
        else; x0 = squeeze(MLE(m,:,iExp-1)); end
        [MLE(m,:,iExp),fMLE(m,:,iExp)] = B{m}.maximizeLikelihood(x0);
    end
end

%% Compute FIM
cellCounts = A.ssaOptions.nSimsPerExpt*ones(size(A.tSpan));  % Number of cells in each experiment.
for m=1:3
    fimResults{m} = B{m}.computeFIM(sensSoln.sens); 
    % Compute the FIM for full observations and no distortion.
    [FIM{m},sFIMcov{m},fimMetrics{m}] = B{m}.evaluateExperiment(fimResults{m},cellCounts);
    fimFreePars = FIM{m}{1}(B{m}.fittingOptions.modelVarsToFit,B{m}.fittingOptions.modelVarsToFit);
    B{m}.makeMleFimPlot(squeeze(MLE(m,:,:)),fimFreePars,[1,2],0.95)
end

%% Optimize Experiment
A.tSpan = [0:1:99]; % Set times at which to compute distributions
A.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = A.solve;  % Solve the sensitivity problem
cellCounts = 100*ones(size(A.tSpan));
nCellsTotal = sum(cellCounts);
% Compute FIM with all parameters free to vary
A.fittingOptions.pdoVarsToFit = 'all';
fimResults = A.computeFIM(sensSoln.sens); 
FIMintuit = A.evaluateExperiment(fimResults,cellCounts);
% Optimize cell counts for different objectives
nCellsOptDet = A.optimizeCellCounts(fimResults,nCellsTotal,'Determinant');
FIMoptDet= A.evaluateExperiment(fimResults,nCellsOptDet);
nCellsOpt1234 = A.optimizeCellCounts(fimResults,nCellsTotal,'[1:4]');
FIMoptDet1234= A.evaluateExperiment(fimResults,nCellsOpt1234);

% Plot results
pars = [[A.parameters{:,2}],[0.7,0.7]]';
for i=1:5
    for j = i+1:6
        subplot(5,5,(i-1)*5+j-1)
        A.makeMleFimPlot([],FIMintuit{1},[j,i],0.95,1,pars); hold on
        A.makeMleFimPlot([],FIMoptDet{1},[j,i],0.95,1,pars)
        A.makeMleFimPlot([],FIMoptDet1234{1},[j,i],0.95,1,pars)
        legend("off")
    end
end
legend('Intuitive','D-Optimality','$\min \det(\Sigma_{\rm model})$','interpreter','latex')
