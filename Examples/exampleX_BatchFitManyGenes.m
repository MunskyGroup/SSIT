addpath(genpath('../src'));

%% Define Base Model Combination
Model_Template = SSIT;
Model_Template.species = {'onGene';'rna'};
Model_Template.initialCondition = [0;0];
Model_Template.propensityFunctions = {'(kon_0+kon_1*Iupstream)*(2-onGene)';...
    'koff_0/(1+akoff)*onGene';'kr_0*(2-onGene)+kr_1*onGene';'gr*rna'};
Model_Template.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model_Template.inputExpressions = {'Iupstream','1+a1*exp(-r1*t*(t>=0))*(1-exp(-r2*t*(t>=0)))'};
Model_Template.parameters = ({'a1',100;'r1',0.01;'r2',0.1;...
    'kon_0',0.01;'kon_1',0.01;'koff_0',20;'akoff',0.2;'kr_0',1;'kr_1',100;'gr',1});
Model_Template.fspOptions.initApproxSS = true;

Model_Template.fittingOptions.modelVarsToFit = 1:10;
Model_Template.fittingOptions.logPrior = @(x)-sum(log10(x).^2/2);
% We generate functions for model propensities
[~,~,Model_Template] = Model_Template.solve;

% Set PDO to Binomial for RNA (assuming 95% drop-out rate).
Model_Template.pdoOptions.type = 'Binomial';
Model_Template.pdoOptions.unobservedSpecies = 'onGene';
Model_Template.pdoOptions.props.CaptureProbabilityS1 = 0;    % Gene State is not measured
Model_Template.pdoOptions.props.CaptureProbabilityS2 = 0.05; % 95% drop out from RNA 
[~,Model_Template] = Model_Template.generatePDO;

%% Load and fit representative data set to get better first parameter guess.
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
Model_Template = Model_Template.loadData(DataFileName,{'rna','DUSP1'});
for i=1:5; [~,~,~,Model_Template] = Model_Template.maximizeLikelihood; end
% Model_Template.makeFitPlot;  % Cluster may crash if no display is set.

%% Generate library of individual gene models
% Specify datafile name and species linking rules
% DataFileName = 'data/ReducedSeqData.csv';
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
TAB = readtable(DataFileName);
geneNames = fields(TAB);
geneNames = geneNames(2:end-4); % The Genes are in Columns 2 -> N-4

for iGene = 1:length(geneNames)
    linkedSpecies = {'rna',geneNames{iGene}};
    Model = Model_Template.loadData(DataFileName,linkedSpecies);
    modelName = ['Model_',geneNames{iGene}];
    eval([modelName,'= Model;']);
    save(['seqModels/',modelName],modelName)
end

%% Fit a multi-model for the 4 genes
% Here we use the multi-model approach on a few genes to constrain the
% parameters of the upstream input signal.

% Select which models to include in multimodel.
Models = {Model_DUSP1,Model_FKBP5,Model_BIRC3,Model_TSC22D3};

% Define how parameters are assigned to sub-models.  All genes are
% assumed to use the same upstream signal, but have different gene bursting
% parameters. 
ParInds = {[1:3,4:10],[1:3,7*1+(4:10)],[1:3,7*2+(4:10)],[1:3,7*3+(4:10)]};

% Define constraint on model parameters. They should be of similar
% magnitudes for each gene unless otherwise demanded by the differences in
% the data.
Constraint = @(x) -var(log10([x(4:10);x(7*1+(4:10));x(7*2+(4:10));x(7*3+(4:10))]));

% Create and initialize multimodel
combinedModel = SSITMultiModel(Models, ParInds, Constraint);
combinedModel = combinedModel.initializeStateSpaces();

%% Fit the multimodel.
% Because there are a lot of parameters, this could take a few rounds to
% get a good MLE. You should be able to get a best posterior of about XXX.
for i = 1:10; [~,~,~,combinedModel] = combinedModel.maximizeLikelihood; end

%% Update upstream signal in all models, fix those parameters, and save.
for iGene = 1:length(geneNames)
    modelName = ['Model_',geneNames{iGene}];
    eval([modelName,'= Model;']);
    eval([modelName,'.parameters(:,2) = num2cell(combinedModel.parameters(1:10));']);
    eval([modelName,'.fittingOptions.modelVarsToFit = 4:10';]);
    save(['seqModels/',modelName],modelName)
end

%% Call Pipeline to Fit Model
% Specify pipeline to apply to model and arguments
% ("../../SSIT/src/exampleData/examplePipelines/fittingPipelineExample.m") 
Pipeline = 'fittingPipelineExample';
pipelineArgs.maxIter = 1000;
pipelineArgs.display = 'iter';
pipelineArgs.makePlot = false;
pipelineArgs.nRounds = 5;

%% Launch cluster jobs for all genes.
for iGene = 1:length(geneNames)
    modelName = ['Model_',geneNames{iGene}];
    saveName = ['seqModels/',modelName];
    % SSIT(saveName,modelName,[],Pipeline,pipelineArgs,saveName);
    
    logfile = ['logFiles/log',modelName];
    % cmd = SSIT.generateCommandLinePipeline(saveName,modelName,[],Pipeline,pipelineArgs,saveName,logfile)
    % cmd = SSIT.generateCommandLinePipeline(saveName,modelName,[],Pipeline,pipelineArgs,saveName,logfile,1)
    cmd = SSIT.generateCommandLinePipeline(saveName,modelName,[],Pipeline,pipelineArgs,saveName,logfile,1,1);
end



