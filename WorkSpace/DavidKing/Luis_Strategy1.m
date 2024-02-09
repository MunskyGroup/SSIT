%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
  addpath(genpath('../../src'));
  DATAFILE = 'datasetDavidKingJan2024';
%% Define SSIT Model
  Model0 = SSIT;
  Model0.species = {'elt2';... % endogenous ELT-2. Unobserved.
                   'elt7';... % endogenous ELT-7. Unobserved.
                   'gfp'};   % GFP, run off of elt-2 promoter. This is our observed variable.
  Model0.initialCondition = [0;0;0];
  Model0.propensityFunctions = {
      'kelt2*(1 + elt7)/(promoter_alpha+elt7)*promoter_beta/(promoter_beta+elt2)'; ... % ELT-2
                                                                'gelt2*elt2';... % degradation
                                                                'kelt7';...      % ELT-7
                                                               'gelt7*elt7';...  % degradation
        'kgfp*(1 + elt7)/(promoter_alpha+elt7)*promoter_beta/(promoter_beta+elt2)';... % GFP                          
                                                               'ggfp*gfp'}; %GFP
  Model0.stoichiometry = [1,-1,0,0,0,0;...
                         0,0,1,-1,0,0;...
                         0,0,0,0,1,-1];

  Model0.parameters = ({ ...
      'kelt2',1;...
      'kelt7',1;...
      'kgfp',1;...
      'gelt2',1;...
      'gelt7',1;...
      'ggfp',16;...
      'promoter_alpha',5; ...
      'promoter_beta',10});

%% set FSP and other settings
disp("Setting parameters.")
Model0.fspOptions.initApproxSS = true;
Model0.tSpan=[0,100];

Model0.solutionScheme = 'FSP';
Model0.fspOptions.fspTol = 1e-3;
disp('Writing propensity functions')
Model0 = Model0.formPropensitiesGeneral('DavidKing');  

[~,Model0.fspOptions.bounds] = Model0.solve;
Model0.fspOptions.bounds(6) = max(30,Model0.fspOptions.bounds(6));

%% Load and Fit GFP reporter data
disp("Loading WT data.");
Model1 = Model0.loadData(DATAFILE,{'gfp','IntensityIntegers'},...
  {'Genotype','ELT-2-GFP';'RNAi','L4440'});
Model1.fittingOptions.modelVarsToFit = [1:8];

disp("Loading elt7 KO data.");
Model2 = Model0.loadData(DATAFILE,{'gfp','IntensityIntegers'; 'elt7', 'zeroes'},...
 {'Genotype','ELT-7-KO';'RNAi','L4440'}); 
Model2.parameters(2,:) = {'kelt7',0};
% define parameters that are in model 2.
M2parInds = [1,3,4,6,8];
Model2.fittingOptions.modelVarsToFit = M2parInds; 

%%
[~,Model1.fspOptions.bounds] = Model1.solve;
[~,Model2.fspOptions.bounds] = Model2.solve;

mu = zeros(size(Model0.parameters,1),1);  % log means
sig2 = ones(size(Model0.parameters,1),1);

% Model1.fittingOptions.logPrior = @(x)sum(-(log10(x)-mu).^2./(2*sig2));
% Model2.fittingOptions.logPrior = @(x)sum(-(log10(x)-mu([1,3,4,6,8])).^2./(2*sig2(M2parInds)));

%% setting up the fit
% initial guess
pars0 = [Model0.parameters{:,2}];
fitOptions = optimset('Display','iter','MaxIter',100);
obj = @(x)-(Model1.computeLikelihood(exp(x))+Model2.computeLikelihood(exp(x(M2parInds))));
initial_error = obj(log(pars0))

%% fit
disp("Fitting...");

% load ParsBrian pars0
for iRound = 1:5
    pars0 = exp(fminsearch(obj,log(pars0),fitOptions))
    
    Model1.parameters(:,2) = num2cell(pars0);
    Model2.parameters(M2parInds,2) = num2cell(pars0(M2parInds));

    [soln1,Model1.fspOptions.bounds] = Model1.solve;
    Model1.fspOptions.bounds(6) = max(30,Model1.fspOptions.bounds(6));
    [~,~,soln1] = Model1.computeLikelihood([],soln1.stateSpace);

    [soln2,Model2.fspOptions.bounds] = Model2.solve;
    [~,~,soln2] = Model2.computeLikelihood([],soln2.stateSpace);
    Model2.fspOptions.bounds(6) = max(30,Model2.fspOptions.bounds(6));

    Model1.makeFitPlot(soln1,1);
    Model2.makeFitPlot(soln2,1);
end

% save ParsBrian pars0



% Model0.parameters(Model0.fittingOptions.modelVarsToFit,2) = num2cell(Model0.maximizeLikelihood([],fitOptions));

%% plot
% disp("Plotting...");
% [fspSoln,Model0.fspOptions.bounds] = Model0.solve;
% Model0.makeFitPlot([],1); % makeSeparatePlotOfData last worked on commit 20a9bbe, not subsequent: 0055cc
% Error with 0055cc:
% Attempt to grow array along ambiguous dimension.
% 
% Error in makeSeparatePlotOfData (line 87)
%                 H1(end+1)=0;

