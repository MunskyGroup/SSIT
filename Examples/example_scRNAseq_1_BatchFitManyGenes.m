% scRNA-seq example 
% * Define and solve template model using FSP
% * Generate Binomial PDO for missing RNA counts
% * Show example for using SSITMultiModel to fit multiple models with some 
%   shared parameters and some distinct parameters
addpath(genpath('../src'));

%% Define Base Model Combination
Model_Template = SSIT;
Model_Template.species = {'onGene';'rna'};
Model_Template.initialCondition = [0;0];
Model_Template.propensityFunctions = ...
    {'(kon_0+kon_1*Iupstream)*(2-onGene)';...
    'koff_0/(1+akoff*Iupstream)*onGene';...
    'kr_0*(2-onGene)+kr_1*onGene';'gr*rna'};
Model_Template.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model_Template.inputExpressions = {'Iupstream',...
                                'exp(-r1*t*(t>=0))*(1-exp(-r2*t*(t>=0)))'};
Model_Template.parameters = ({'r1',0.01; 'r2',0.1; 'kon_0',0.01;...
                              'kon_1',0.01; 'koff_0',20; 'akoff',0.2;...
                              'kr_0',1; 'kr_1',100; 'gr',1});

Model_Template.fspOptions.initApproxSS = true;
Model_Template.fittingOptions.modelVarsToFit = 1:9;
Model_Template.fittingOptions.logPrior = @(x)-sum(log10(x).^2/2);

% We generate functions for model propensities
Model_Template = Model_Template.formPropensitiesGeneral('Model_Template');
[~,~,Model_Template] = Model_Template.solve;

%% Set PDO to Binomial for RNA (assuming 95% drop-out rate)
Model_Template.pdoOptions.type = 'Binomial';
Model_Template.pdoOptions.unobservedSpecies = 'onGene';
Model_Template.pdoOptions.props.CaptureProbabilityS1 = 0;    % Gene State is not measured
Model_Template.pdoOptions.props.CaptureProbabilityS2 = 0.05; % 95% drop out from RNA
[~,Model_Template] = Model_Template.generatePDO();

%% Load and fit representative data set to get better first parameter guess
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
Model_Template = Model_Template.loadData(DataFileName,{'rna','DUSP1'});
for i=1:5
    [~,~,~,Model_Template] = Model_Template.maximizeLikelihood;
end
Model_Template.makeFitPlot;  % Cluster may crash if no display is set.

%% Generate library of individual gene models
% Specify datafile name and species linking rules
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
TAB = readtable(DataFileName);
geneNames = fields(TAB);
geneNames = geneNames(2:end-4); % The Genes are in Columns 2 -> N-4

if ~exist('seqModels','dir'); mkdir('seqModels'); end
for iGene = 1:length(geneNames)
    linkedSpecies = {'rna',geneNames{iGene}};
    Model = Model_Template.loadData(DataFileName,linkedSpecies);
    modelName = ['Model_',geneNames{iGene}];
    assignin('base',modelName,Model);
    save(['seqModels/',modelName],modelName);
end

%% Fit a multi-model for the 4 genes
% Here we use the multi-model approach on a few genes to constrain the
% parameters of the upstream input signal.

% Select which models to include in multimodel.
Models = {Model_DUSP1,Model_RUNX1,Model_BIRC3,Model_TSC22D3};
modelNames = {'Model_DUSP1','Model_RUNX1','Model_BIRC3','Model_TSC22D3'};

% Define how parameters are assigned to sub-models.  All genes are
% assumed to use the same upstream signal, but have different gene bursting
% parameters.
ParInds = {[1,2,3:9],[1,2,7*1+(3:9)],[1,2,7*2+(3:9)],[1,2,7*3+(3:9)]};

% Define constraint on model parameters. They should be of similar
% magnitudes for each gene unless otherwise demanded by the differences in
% the data.
Constraint = @(x) -var(log10([x(3:9);x(7*1+(3:9));x(7*2+(3:9));x(7*3+(3:9))]));

% Create and initialize multimodel
combinedModel = SSITMultiModel(Models, ParInds, Constraint);
combinedModel = combinedModel.initializeStateSpaces();

%% Fit the multimodel.
% Because there are a lot of parameters, this could take a few rounds to
% get a good MLE. You should be able to get a best posterior of about XXX.
fitOptions = optimset('Display','iter','MaxIter',1000);
fitOptions.suppressExpansion = true;
for i = 1:10
    [~,~,~,combinedModel] = combinedModel.maximizeLikelihood([],fitOptions);
    save('seqModels/CombinedModel4Genes','combinedModel');
end

%% Update upstream signal in all models, fix those parameters, and save
for iGene = 1:length(geneNames)
    modelName = ['Model_',geneNames{iGene}];
    eval([modelName,'.parameters(:,2) = num2cell(combinedModel.parameters(1:9));']);
    eval([modelName,'.fittingOptions.modelVarsToFit = 3:9;';]);
    save(['seqModels/',modelName],modelName);
end

%% Update the four models used in multimodel demo with multimodel-fitted
%  parameters (Model_DUSP1, Model_RUNX1, Model_BIRC3, Model_TSC22D3):
for iGene = 1:4
    modelName = modelNames{iGene};

    Models{iGene}.parameters(:,2) = ...
        combinedModel.SSITModels{iGene}.parameters(:,2);

    % Update the workspace variable Model_XXXX
    assignin('base', modelName, Models{iGene});

    % Save Model_XXXX into seqModels/Model_XXXX.mat
    m = struct();
    m.(modelName) = Models{iGene};
    save(fullfile('seqModels',[modelName '.mat']),'-struct','m',modelName);

    % Make fit plots:
    Models{iGene}.makeFitPlot
end