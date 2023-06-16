%% example_BPPB -- example script to show SSIT for the BPPB Seminar
clear all
close all
clc
addpath('../')

%% test 2
Model = SSIT;
Model.species = {'x1';'x2'};
Model.initialCondition = [0;0];
Model.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
Model.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model.inputExpressions = {'IGR','1+(t>=0)*a1*exp(-r1*t)*(1-exp(-r2*t))'};
Model.parameters = ({'koff',0.14;'kon',0.14;'kr',5;'gr',0.05;...
                  'a1',0.4;'r1',0.04;'r2',0.1});
Model.fspOptions.initApproxSS = true;
Model.fittingOptions.modelVarsToFit = ones(1,7,'logical');
Model.tSpan = [0 10 20 30 40 50 60 75 90 120 150 180]; 
 
simple_Model = load('simple_dusp1_model.mat').simple_Model;
Model.parameters=simple_Model.parameters;

[FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis
[FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis

%%
Model.ssaOptions.nSimsPerExpt = 100;
Model.ssaOptions.Nexp = 50; 
Model.fspOptions.fspTol = 1e-8;  % Set FSP error tolerance.
Model.sampleDataFromFSP(FSPsoln,'simple_dusp1_model_testC.csv'); 

%% Solve Sensitivity using FSP
Model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
Model.sensOptions.solutionMethod = 'finiteDiff';
[sensSoln] = Model.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

%% Associate Datasets with FSP Models.
Model.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
Model = Model.loadData('simple_dusp1_model_testC.csv',{'x2','exp1_s2'});
Model.initialTime = 0;
Model.fittingOptions.timesToFit = ones(1,length(Model.tSpan),'logical');
Model.makeFitPlot

%% Find MLE for each simulated data set.
Model.fittingOptions.modelVarsToFit = [1:4];

par2Fit=length(find(Model.fittingOptions.modelVarsToFit));
fitOptions = optimset('Display','iter','MaxIter',200);
for i=1:3
    if i==1
        MLE = zeros(1,par2Fit,Model.ssaOptions.Nexp); fMLE = inf(1,par2Fit,Model.ssaOptions.Nexp);
    end
    parfor iExp = 1:Model.ssaOptions.Nexp
        tmpModel = Model.loadData('simple_dusp1_model_testC.csv',{'x1',['exp',num2str(iExp),'_s1'];'x2',['exp',num2str(iExp),'_s2']}); % Link non-distorted data.
        %        tmpModel = Model.loadData('simple_dusp1_model_testC.csv',{'x2',['exp',num2str(iExp),'_s2']}); % Link non-distorted data.

        % Link non-distorted data.
        if i==1
            x0 = [tmpModel.parameters{tmpModel.fittingOptions.modelVarsToFit,2}]';
        else
            x0 = squeeze(MLE(1,:,iExp));
        end
        [MLE(1,:,iExp),fMLE(1,:,iExp)] = tmpModel.maximizeLikelihood(x0,fitOptions);
    end

    %% Compute FIM
    cellCounts = Model.ssaOptions.nSimsPerExpt*ones(size(Model.tSpan));  % Number of cells in each experiment.
    fimResults = Model.computeFIM(sensSoln.sens);
    Model.pdoOptions.unobservedSpecies = [];
    % Compute the FIM for full observations and no distortion.
    [FIM,sFIMcov,fimMetrics] = Model.evaluateExperiment(fimResults,cellCounts);
    fimFreePars = FIM(Model.fittingOptions.modelVarsToFit,Model.fittingOptions.modelVarsToFit);
    figure(1);clf;
    for iy = 1:par2Fit-1
        for ix = iy+1:par2Fit
            subplot(par2Fit-1,par2Fit-1,(par2Fit-1)*(iy-1)+ix-1);
            Model.makeMleFimPlot(squeeze(MLE(1,:,:)),fimFreePars,[ix,iy],0.95)
        end
    end
    drawnow
end

%% Optimize Experiment
% A.tSpan = [0:1:99]; % Set times at which to compute distributions
% A.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
% [sensSoln,bounds] = A.solve;  % Solve the sensitivity problem
% cellCounts = 100*ones(size(A.tSpan));
% nCellsTotal = sum(cellCounts);
% % Compute FIM with all parameters free to vary
% A.fittingOptions.pdoVarsToFit = 'all';
% fimResults = A.computeFIM(sensSoln.sens); 
% FIMintuit = A.evaluateExperiment(fimResults,cellCounts);
% % Optimize cell counts for different objectives
% nCellsOptDet = A.optimizeCellCounts(fimResults,nCellsTotal,'Determinant');
% FIMoptDet= A.evaluateExperiment(fimResults,nCellsOptDet);
% nCellsOpt1234 = A.optimizeCellCounts(fimResults,nCellsTotal,'[1:4]');
% FIMoptDet1234= A.evaluateExperiment(fimResults,nCellsOpt1234);
% 
% % Plot results
% pars = [[A.parameters{:,2}],[0.7,0.7]]';
% for i=1:5; for j = i+1:6
%         subplot(5,5,(i-1)*5+j-1)
%         A.makeMleFimPlot([],FIMintuit,[j,i],0.95,1,pars([j,i])); hold on
%         A.makeMleFimPlot([],FIMoptDet,[j,i],0.95,1,pars([j,i]))
%         A.makeMleFimPlot([],FIMoptDet1234,[j,i],0.95,1,pars([j,i]))
%         legend("off")
% end;end
% legend('Intuitive','D-Optimality','$\min \det(\Sigma_{\rm model})$','interpreter','latex')
