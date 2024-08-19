%% example_DUSP1_Regulation
% Example script to show SSIT application to a simple bursting gene
% expression model that ios being fit to smFISH data for DUSP1 upon
% stimulation with Dexamethsone.
clear all
clc
close all
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
Model1.makePlot(Mod1SensSoln,'marginals',[],false,[4,5,6],{'linewidth',2}) % Plot marginal sensitivities

%% STEP4 == Define Binomial Probabilistic Distortion Operator
Model2 = Model1;  % Make a copy of the original model
Model2.pdoOptions.type = 'Binomial';
Model2.pdoOptions.props.CaptureProbabilityS1 = 0;  % Distortion for OFF species (unobserved)
Model2.pdoOptions.props.CaptureProbabilityS2 = 0;  % Distortion for ON species (unobserved)
Model2.pdoOptions.props.CaptureProbabilityS3 = 0.7;% Distortion for RNA species
Model2.pdoOptions.PDO = Model2.generatePDO(Model2.pdoOptions,[0,0,0.7],Mod1FSPsoln.fsp,true);
figure(20); contourf(log10(Model2.pdoOptions.PDO.conditionalPmfs{3})); colorbar
xlabel('"true" number of mRNA'); ylabel('observed number of mRNA'); set(gca,'fontsize',15);

%% STEP5 == Apply PDO to FSP and Sensitivity Calculations
Model2.solutionScheme = 'FSP'; % Set solution scheme to FSP.
Model2.makePlot(Mod1FSPsoln,'marginals',[],true,[1,2,3],{'linewidth',2})  % Plot Distorted Marginals
Model2.solutionScheme = 'fspSens'; % Set solution scheme to Sensitivity
Model2.makePlot(Mod1SensSoln,'marginals',[],true,[4,5,6],{'linewidth',2})    % Plot Distorted Sensitivities

%% STEP6 == Generate Simulated Data for Fitting
Model2.ssaOptions.Nexp = 10;   % Number of independent data sets to generate.
Model2.ssaOptions.nSimsPerExpt = 200; % Number of cells to include at each time point for each data set.
Model2.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
Model2.sampleDataFromFSP(Mod1FSPsoln,'DUSP1SSAData50Expts.csv')  % Generate and save all data.

%% STEP7 == Associate Datasets with FSP Models
Model2.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
B{1} = Model2.loadData('DUSP1SSAData50Expts.csv',{'rna','exp1_s3'});
B{1}.pdoOptions.PDO = [];  % Do not use PDO.
B{2} = Model2.loadData('DUSP1SSAData50Expts.csv',{'rna','exp1_s3_Distorted'});
B{2}.pdoOptions.PDO = [];  % Do not use PDO.
B{3} = Model2.loadData('DUSP1SSAData50Expts.csv',{'rna','exp1_s3_Distorted'});
% B{3}.pdoOptions.unobservedSpecies = {'offGene';'onGene'};

%% STEP8 == Sweep overParameters and Plot Likelihood Functions
fitErrorsB1 = B{1}.likelihoodSweep([3,4],linspace(.5,1.5,11),true);
title('Ideal data. Original FSP.')
fitErrorsB2 = B{2}.likelihoodSweep([3,4],linspace(.5,1.5,11),true);
title('Binomial data distortion. Original FSP.')
fitErrorsB3 = B{3}.likelihoodSweep([3,4],linspace(.5,1.5,11),true);
title('Binomial data distortion. FSP+PDO.')

%% STEP9 == Find MLEs for Simulated Datasets
clear MLE fMLE
for iExp = 1:Model2.ssaOptions.Nexp
    for m = 1:3
        switch m
            case 1
                b = B{1}.loadData('DUSP1SSAData50Expts.csv',{'rna',['exp',num2str(iExp),'_s3']}); % Link non-distorted data.
            case 2
                b = B{2}.loadData('DUSP1SSAData50Expts.csv',{'rna',['exp',num2str(iExp),'_s3_Distorted']}); % Link non-distorted data.
            case 3
                b = B{3}.loadData('DUSP1SSAData50Expts.csv',{'rna',['exp',num2str(iExp),'_s3_Distorted']}); % Link non-distorted data.
        end
        b.fittingOptions.modelVarsToFit = [2,3];  
        x0 = [b.parameters{b.fittingOptions.modelVarsToFit,2}]';
        [MLE(m,:,iExp),fMLE(m,iExp)] = b.maximizeLikelihood(x0);
        [iExp;m;squeeze(MLE(m,:,iExp))']
    end

end

%% STEP10 == Compute FIM and Compare to MLE Spread
for m=1:3
    fimResults{m} = B{m}.computeFIM; 
    % Compute the FIM for full observations and no distortion.
    FIM{m} = B{m}.evaluateExperiment(fimResults{m},B{1}.dataSet.nCells);
    fimFreePars = FIM{m}{1}(b.fittingOptions.modelVarsToFit,b.fittingOptions.modelVarsToFit);
    b.makeMleFimPlot(squeeze(MLE_tmp(m,:,:)),fimFreePars,[1,2],0.95,100+m,[0.002,1])
end

%% STEP11 == Optimize Experiment (no Distortion)
Model3=Model1;            % Make copy of original model
Model3.tSpan = [0:1:180]; % Set allowable experiment times
nCellsIntuitive = 100*ones(size(Model3.tSpan)); % Set Intuitive Expt Design
nCellsTotal = sum(nCellsIntuitive); % Compute total number of cells

Model3.fittingOptions.modelVarsToFit = 'all'; % Allow Model parameters to be free
fimResults = Model3.computeFIM; % Compute FIM at all Times

% Optimize cell counts for different objectives
nCellsOptDet  = Model3.optimizeCellCounts(fimResults,nCellsTotal,'Determinant');
nCellsOpt1234 = Model3.optimizeCellCounts(fimResults,nCellsTotal,'[1:4]');

% Compute total FIM under different designs
FIMintuit = Model3.evaluateExperiment(fimResults,nCellsIntuitive);
FIMoptDet = Model3.evaluateExperiment(fimResults,nCellsOptDet);
FIMoptDet1234= Model3.evaluateExperiment(fimResults,nCellsOpt1234);

% Plot results
pars = Model3.parameters{:,2};
for i=1:6
    for j = i+1:7
        subplot(6,6,(i-1)*6+j-1)
        Model3.makeMleFimPlot([],FIMintuit{1},[j,i],0.95,1,pars); hold on
        Model3.makeMleFimPlot([],FIMoptDet{1},[j,i],0.95,1,pars)
        Model3.makeMleFimPlot([],FIMoptDet1234{1},[j,i],0.95,1,pars)
        legend("off")
    end
end
legend('Intuitive','D-Optimality','$\min \det(\Sigma_{\rm DUSP1})$','interpreter','latex')
%% STEP12 == Optimize Experiment (including Distortion)
Model4=Model2;            % Make copy of original model
Model4.tSpan = [0:1:180]; % Set allowable experiment times
nCellsIntuitive = 100*ones(size(Model4.tSpan)); % Set Intuitive Expt Design
nCellsTotal = sum(nCellsIntuitive); % Compute total number of cells

Model4.fittingOptions.modelVarsToFit = 'all'; % Allow Model parameters to be free
Model4.fittingOptions.pdoVarsToFit = [3]; % Allow PDO parameter (RNA) to be free
fimResults = Model4.computeFIM; % Compute FIM at all Times

% Optimize cell counts for different objectives
nCellsOptDet  = Model4.optimizeCellCounts(fimResults,nCellsTotal,'tr[1:7,10]');

% Compute total FIM under different designs
FIMintuit = Model4.evaluateExperiment(fimResults,nCellsIntuitive);
FIMoptDet = Model4.evaluateExperiment(fimResults,nCellsOptDet);

% Plot results
pars = [[Model4.parameters{:,2}],0.7];
for i=1:7
    for j = i+1:8
        subplot(7,7,(i-1)*7+j-1)
        F = FIMintuit{1}([1:7,10],[1:7,10]);
        Model4.makeMleFimPlot([],F,[j,i],0.95,1,pars); hold on

        F = FIMoptDet{1}([1:7,10],[1:7,10]);
        Model4.makeMleFimPlot([],F,[j,i],0.95,1,pars)

        legend("off")
    end
end
legend('Intuitive','D-Optimality','$\min \det(\Sigma_{\rm DUSP1})$','interpreter','latex')
