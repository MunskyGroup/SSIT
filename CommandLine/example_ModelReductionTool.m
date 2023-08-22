clear all
close all
clc

testModel = 3;
% redType = {'POD 2nd',100};
% redType = {'Log Lump QSSA',40};
redType = {'QSSA', 1};
% redType = {'No Transform',40};
%  redType = {'Proper Orthogonal Decomposition',20};
% redType = {'Eigen Decomposition Initial',30};
% Define SSIT Model
% SSIT models are defined as usual:
switch testModel
    case 1 % Poisson Process
        Model1 = SSIT;
        Model1.species = {'x1'};
        Model1.initialCondition = [0];
        Model1.propensityFunctions = {'kr';'gr*x1'};
        Model1.stoichiometry = [1,-1];
        Model1.parameters = ({'kr',40;'gr',1});
        Model1.fspOptions.initApproxSS = false;
        Model1.tSpan = linspace(0,5,10);
    case 2 % Poisson Start at SS.
        Model1 = SSIT;
        Model1.species = {'x1'};
        Model1.initialCondition = [0];
        Model1.propensityFunctions = {'kr';'gr*x1'};
        Model1.stoichiometry = [1,-1];
        Model1.parameters = ({'kr',40;'gr',1});
        Model1.fspOptions.initApproxSS = true;
        Model1.tSpan = linspace(0,5,10);
    case 3 % Two Species Poisson.
        Model1 = SSIT;
        Model1.species = {'x1';'x2'};
        Model1.initialCondition = [0;0];
        Model1.propensityFunctions = {'kr';'gr*x1';'kp';'gp*x2'};
        Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
        Model1.parameters = ({'kr',40;'gr',1;'kp',20;'gp',1});
        Model1.fspOptions.initApproxSS = false;
        Model1.tSpan = linspace(0,5,12);
    
    case 4 % Time varying model (DUSP1)
        Model1 = SSIT;
        Model1.species = {'ActiveGene';'mRNA'};
        Model1.initialCondition = [0;0];
%         Model1.propensityFunctions = {'kon*(2-x1)';'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
%         Model1.stoichiometry = [1,1,-1,0,0;0,0,0,1,-1];
%         Model1.inputExpressions = {'IGR','a1*exp(-r1*t)*(1-exp(-r2*t))'};

        Model1.propensityFunctions = {'kon*(1+IGR)*(2-ActiveGene)';'koff*ActiveGene';'kr*ActiveGene';'gr*mRNA'};
        Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
        Model1.inputExpressions = {'IGR','a1*exp(-r1*t)*(1-exp(-r2*t))'};

        %         Model1.propensityFunctions = {'kon*(2-x1)';'kon*IGR*(2-x1)';'-kon*IGR*x1';'koff*x1';'kr*x1';'gr*x2'};
        %         Model1.stoichiometry = [1,1,1,-1,0,0;0,0,0,0,1,-1];
        Model1.parameters = ({'koff',0.14;'kon',0.14;'kr',10;'gr',0.01;...
            'a1',0.4;'r1',0.04;'r2',0.1});
        Model1.fspOptions.initApproxSS = true;
        %   Model1.tSpan = linspace(0,180,20);
        Model1.tSpan = linspace(0,180,12);
%         Model1.fspOptions.bounds(4) = 100;
        Model1.fspOptions.verbose = true;
    case 5 % Time varying Poisson Process
        Model1 = SSIT;
        Model1.species = {'x1'};
        Model1.initialCondition = [0];
        Model1.propensityFunctions = {'kr+I1';'gr*x1'};
        Model1.stoichiometry = [1,-1];
        Model1.parameters = ({'kr',40;'gr',1});
        Model1.inputExpressions = {'I1','40*t'};
        Model1.fspOptions.initApproxSS = false;
        Model1.tSpan = linspace(0,5,10);
end
      


%%
[fspSoln,Model1.fspOptions.bounds] = Model1.solve;
tic
[fspSoln,Model1.fspOptions.bounds] = Model1.solve(fspSoln.stateSpace);
fullModelSolveTime = toc
Model1.makePlot(fspSoln,'meansAndDevs',[],[],1)
Model1.makePlot(fspSoln,'marginals',[],[],[2,3])
% Model1.makePlot(fspSoln,'joints',[1:100:1000],[],[12])
%%
for iMR = 1:size(redType,1)
    Model2{iMR} = Model1;
    Model2{iMR}.fspOptions.fspTol = inf;
    Model2{iMR}.modelReductionOptions.reductionType = redType{iMR,1};
    Model2{iMR}.modelReductionOptions.reductionOrder = redType{iMR,2};

    if strcmp(redType{1},'QSSA')
        Model2{iMR}.modelReductionOptions.qssaSpecies = redType{iMR,2};
    end

    Model2{iMR} = Model2{iMR}.computeModelReductionTransformMatrices(fspSoln);

    Model2{iMR}.modelReductionOptions.useModReduction = true;
    
    tic
    fspSoln2{iMR} = Model2{iMR}.solve(fspSoln.stateSpace);
    redModelSolveTime = toc

    Model2{iMR}.makePlot(fspSoln2{iMR},'meansAndDevs',[],[],1)
    Model2{iMR}.makePlot(fspSoln2{iMR},'marginals',[],[],[2,3])
%     Model2{iMR}.makePlot(fspSoln2{iMR},'joints',[1:100:1000],[],[13])
end
return
%%
Model3 = Model2{1};
Model3.modelReductionOptions.reductionType ='POD Update';
Model3 = Model3.computeModelReductionTransformMatrices(fspSoln);


%% Fitting using reduced model

Model2{1}.tSpan = 0;
Model2{1} = Model2{1}.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});
Model2{1}.fspOptions.fspTol = inf;
Model2{1}.fittingOptions.modelVarsToFit = 1:2;
fitOptions = optimset('Display','iter','MaxIter',4000);
%
Model2{1}.parameters(Model2{1}.fittingOptions.modelVarsToFit,2) = num2cell(Model2{1}.maximizeLikelihood([],fitOptions));
Model2{1}.makeFitPlot;
