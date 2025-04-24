clear
addpath(genpath('../../src'))
%% Create GR Model Template
ModelGR = SSIT;
ModelGR.species = {'cytGR';
    'nucGR'};
ModelGR.initialCondition = [20;1];

ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';
    'knc*nucGR';...
    'kg1';
    'gg1*cytGR';
    'gg2*nucGR'};
ModelGR.stoichiometry = [-1,1,1,-1,0;...
    1,-1,0,0,-1];

ModelGR.parameters = ({'koff',0.1; ...
    'kon',0.1; ...
    'kr',1; ...
    'gr',0.02;...
    'kcn0',0.0054; ...
    'kcn1',0.11; ...
    'gDex',4e-6; ...
    'knc',0.017; ...
    'kg1',0.011;...
    'gg1',3.3e-5; ...
    'gg2',0.0041; ...
    'MDex',15; ...
    'Dex0',100});

ModelGR.fspOptions.initApproxSS = true;

ModelGR.fittingOptions.modelVarsToFit = (5:12);

ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
ModelGR = ModelGR.formPropensitiesGeneral('ModelGRonly');
ModelGR.customConstraintFuns = {'cytGR+nucGR'};
save('savedModelGR','ModelGR')

%% Load Data into Model
dataFile = 'data/dataSet.csv';
speciesMapping = {'nucGR','normgrnuc';'cytGR','normgrcyt'};
dataSets = {dataFile,speciesMapping,{'Dex_Conc','1'}};
ModelGRwithData = SSIT('savedModelGR.mat','ModelGR',dataSets);

%% Load Data into Model and Fit
dataFile = 'data/dataSet.csv';
speciesMapping = {'nucGR','normgrnuc';'cytGR','normgrcyt'};
dataSets = {dataFile,speciesMapping,{'Dex_Conc','1'}};
Pipeline = 'fittingPipelineExample';
pipelineArgs.maxIter = 200;
pipelineArgs.display = 'iter';
pipelineArgs.makePlot = false;
saveFile = 'exampleResults.mat';

% Create model from preset, associate with data, run
% 'fittingPipeline', and save result.
ModelGRwithData = SSIT('savedModelGR.mat','ModelGR',dataSets,Pipeline,pipelineArgs,saveFile);

%%
% Load model from file, run 'fittingPipeline', and save result.
pipelineArgs.maxIter = 100;
pipelineArgs.display = 'iter';
pipelineArgs.makePlot = true;
Model = SSIT(saveFile,'ModelGR',[],Pipeline,pipelineArgs,saveFile);
%% Create GR MultiModel
% Load Saved Model from GR Model Template.
ModelGR = SSIT('savedModelGR.mat','ModelGR');

% Make copies for different conditions.
ModelGR1 = ModelGR;
ModelGR1.parameters(13,:)  = {'Dex0', 1};

ModelGR10 = ModelGR;
ModelGR10.parameters(13,:) = {'Dex0', 10};

ModelGR100 = ModelGR;
ModelGR10.parameters(13,:) = {'Dex0', 100};

% Set parameter use (all use same parameters)
MultiModelGRparameterMap = {(1:8),(1:8),(1:8)};

% Set log-prior on parameters
log10PriorMean = [-1 -1 0 -2,...
    -1 -3 -2 -1 -2 -2 -2 0.5, 2];
log10PriorStd = 2*ones(1,13);
logPriorGR = @(x)-sum((log10(x(:))-log10PriorMean(5:12)).^2./(2*log10PriorStd(5:12).^2),"all");

% Make MultiModel
combinedGRModel = SSITMultiModel({ModelGR1,ModelGR10,ModelGR100}, ...
    MultiModelGRparameterMap, ...
    logPriorGR);

% Save MultiModel
save('savedGRmultiModel',"combinedGRModel")

%% Load, Assign data, Fit, and Save GR multimodel
% clear
% Define data sets to use
dataFile = 'data/dataSet.csv';
speciesMapping = {'nucGR','normgrnuc';'cytGR','normgrcyt'};
dataSets = {dataFile,speciesMapping,{'Dex_Conc','1'};...
            dataFile,speciesMapping,{'Dex_Conc','10'};...
            dataFile,speciesMapping,{'Dex_Conc','100'}};

% specify fitting routine
pipeline = 'MultiModelFittingPipeline';
pipelineArgs.runFit = true;
pipelineArgs.fitOpts.maxIter = 20;
pipelineArgs.fitOpts.display = 'iter';
pipelineArgs.makePlots = false;
pipelineArgs.runMH = false;

% Specify save name
saveName = 'GRModelFitResults';

SSIT('savedGRmultiModel','combinedGRModel',dataSets,pipeline,pipelineArgs,saveName);

%% Run Metropolis Hastings Pipeline on MultiModel
pipeline = 'MultiModelFittingPipeline';
pipelineArgs.runFit = false;
pipelineArgs.runMH = true;
pipelineArgs.MHopts = struct('numberOfSamples',100);
pipelineArgs.makePlots = true;

saveName = 'GRModelFitResults';
SSIT('GRModelFitResults','combinedGRModel',{},pipeline,pipelineArgs,saveName);