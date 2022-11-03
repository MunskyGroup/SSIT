clear all
%% Define SSIT Model
    Model = SSIT;
    Model.species = {'x1';'x2'};
    Model.initialCondition = [0;0];
    Model.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
    Model.stoichiometry = [1,-1,0,0;0,0,1,-1];
    Model.inputExpressions = {'IGR','a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'};
    Model.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01; ...
                   'a0',0.006;'a1',0.4;'r1',0.04;'r2',0.1});
    Model.initialTime = -120;  % large negative time to simulate steady state at t=0

%% Load and Fit smFISH Data
    Model = Model.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x2','RNA_nuc'});
    Model.tSpan = unique([Model.initialTime,Model.dataSet.times]);
    fitOptions = optimset('Display','iter','MaxIter',100);
    [pars,likelihood] = Model.maximizeLikelihood([],fitOptions);

%% Update Model and Make Plots of Results
    Model.parameters(:,2) = num2cell(pars);
    Model.makeFitPlot

%% Metropolis Hastings
    MHOptions = struct('numberOfSamples',1000,'burnin',500);
    [pars,likelihood,chainResults] = Model.maximizeLikelihood([],MHOptions,'MetropolisHastings');






