% In this example, we show how to create reduced FSP models using different
% types of projectionbased transformations.
close all 
clear all
addpath('../CommandLine');

%% First, choose a model on which to illustrate the reduction approximation,
% or you can create your own. See below for the codes to create each model.
testModel = 3; 

%% Next, choose which type of model reduction to apply. Options include:
%   'Proper Orthogonal Decomposition' - solve the FSP once and then uses
%       POD to construct a reduced basis set that covers the current FSP
%       solution.  For best use, this reduction should be found using a
%       fine time resolution in the calculation of the FSP. The size of the
%       reduced model must be specified as 'reductionOrder'. Because the
%       POD requires a solution of the FSP, this reduction is usualy only
%       helpful for situations where many solutions are needed (e.g.,
%       during model fitting).
%   'Log Lump QSSA' - forms a coarse rectangular mesh with grid points chosen
%       logarithmically using the current FSP bounds. The number of grid
%       lines must be specified in 'reductionOrder'.
%   'Eigen Decomposition Initial' - reduction to consider only the space
%       spanned by the initial condition plus the eigvenvectors corresponding
%       to the eigenvalues wite the largest real values.  The number of modes
%       to consider in the reduction is specified in 'reductionOrder'. Fort
%       time varying systems, the basis vectors are found using the
%       infinitesimal generator at t=0.  
%   'No Transform' - test case where no reduction is applied.
%   'QSSA' - Reduction using QSSA applied to a specific species or set of
%       species. The list of species to be assumed at QSSA must be
%       specified in a vector 'reductionSpecies'. 
%%
reductionType = 'Proper Orthogonal Decomposition'; %{'Log Lump QSSA','Proper Orthogonal Decomposition','QSSA'};
reductionOrder = 20;
qssaSpecies = 2;
podTimeSetSize = 100;

% Define SSIT Model
% SSIT models are defined as usual:
switch testModel
    case 1 % Poisson Process
        Model1 = SSIT;
        Model1.species = {'x1'};
        Model1.initialCondition = 0;
        Model1.propensityFunctions = {'kr';'gr*x1'};
        Model1.stoichiometry = [1,-1];
        Model1.parameters = ({'kr',100;'gr',1});
        Model1.fspOptions.initApproxSS = false;
        Model1.tSpan = linspace(0,3,10);
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
        Model1.propensityFunctions = {'kon*(1+IGR)*(2-ActiveGene)';'koff*ActiveGene';'kr*ActiveGene';'gr*mRNA'};
        Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
        Model1.inputExpressions = {'IGR','a1*exp(-r1*t)*(1-exp(-r2*t))'};
        Model1.parameters = ({'koff',0.14;'kon',0.14;'kr',10;'gr',0.01;...
            'a1',0.4;'r1',0.04;'r2',0.1});
        Model1.fspOptions.initApproxSS = true;
        Model1.tSpan = linspace(0,180,12);
end
 
[fspSoln,Model1.fspOptions.bounds] = Model1.solve;

if strcmp(reductionType,'Proper Orthogonal Decomposition')
    tSpan = Model1.tSpan;
    Model1.tSpan = linspace(min(Model1.tSpan),max(Model1.tSpan),podTimeSetSize);
    [fspSoln,Model1.fspOptions.bounds] = Model1.solve(fspSoln.stateSpace);
else
    tic
    [fspSoln,Model1.fspOptions.bounds] = Model1.solve(fspSoln.stateSpace);
    fullModelSolveTime = toc
end

if strcmp(reductionType,'Proper Orthogonal Decomposition')
    Model1.tSpan = tSpan;
end

Model2 = Model1;
Model2.modelReductionOptions.useModReduction = true;
Model2.fspOptions.fspTol = inf;
Model2.modelReductionOptions.reductionType = reductionType;
Model2.modelReductionOptions.reductionOrder = reductionOrder;
Model2.modelReductionOptions.qssaSpecies = qssaSpecies;
Model2 = Model2.computeModelReductionTransformMatrices(fspSoln);

tic
fspSoln2 = Model2.solve(fspSoln.stateSpace);
redModelSolveTime = toc

if strcmp(reductionType,'Proper Orthogonal Decomposition')
    tic
    [fspSoln,Model1.fspOptions.bounds] = Model1.solve(fspSoln.stateSpace);
    fullModelSolveTime = toc
end


Model1.makePlot(fspSoln,'meansAndDevs',[],[],1)
Model1.makePlot(fspSoln,'marginals',[],[],[2,3])

Model2.makePlot(fspSoln2,'meansAndDevs',[],[],1)
Model2.makePlot(fspSoln2,'marginals',[],[],[2,3])
return
%%
Model3 = Model2;
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
