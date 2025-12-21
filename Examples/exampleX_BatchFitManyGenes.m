%Multi-Model Example.
addpath(genpath('../src'));

%% Define Base Model Combination
Model_Template = SSIT;
Model_Template.species = {'onGene';'rna'};
Model_Template.initialCondition = [0;0];
Model_Template.propensityFunctions = {'(kon_0+kon_1*Iupstream)*(2-onGene)';...
    'koff_0/(1+akoff*Iupstream)*onGene';'kr_0*(2-onGene)+kr_1*onGene';'gr*rna'};
Model_Template.stoichiometry = [1,-1,0,0;0,0,1,-1];
Model_Template.inputExpressions = {'Iupstream','exp(-r1*t*(t>=0))*(1-exp(-r2*t*(t>=0)))'};
Model_Template.parameters = ({'r1',0.01;'r2',0.1;...
    'kon_0',0.01;'kon_1',0.01;'koff_0',20;'akoff',0.2;'kr_0',1;'kr_1',100;'gr',1});
Model_Template.fspOptions.initApproxSS = true;

Model_Template.fittingOptions.modelVarsToFit = 1:9;
Model_Template.fittingOptions.logPrior = @(x)-sum(log10(x).^2/2);
% We generate functions for model propensities
[~,~,Model_Template] = Model_Template.solve;

% Set PDO to Binomial for RNA (assuming 95% drop-out rate).
% TODO - Alex, please update to your newer version.
Model_Template.pdoOptions.type = 'Binomial';
Model_Template.pdoOptions.unobservedSpecies = 'onGene';
Model_Template.pdoOptions.props.CaptureProbabilityS1 = 0;    % Gene State is not measured
Model_Template.pdoOptions.props.CaptureProbabilityS2 = 0.05; % 95% drop out from RNA 
[~,Model_Template] = Model_Template.generatePDO;

%% Load and fit representative data set to get better first parameter guess.
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
Model_Template = Model_Template.loadData(DataFileName,{'rna','DUSP1'});
for i=1:5
    [~,~,~,Model_Template] = Model_Template.maximizeLikelihood; 
end
% Model_Template.makeFitPlot;  % Cluster may crash if no display is set.

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

%% Update upstream signal in all models, fix those parameters, and save.
for iGene = 1:length(geneNames)
    modelName = ['Model_',geneNames{iGene}];
    eval([modelName,'.parameters(:,2) = num2cell(combinedModel.parameters(1:9));']);
    eval([modelName,'.fittingOptions.modelVarsToFit = 3:9;';]);
    save(['seqModels/',modelName],modelName);
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
    logfile = ['logFiles/log',modelName];
    runNow = true;
    runCluster = true;
    cmd = SSIT.generateCommandLinePipeline(saveName,modelName,[],Pipeline,...
        pipelineArgs,saveName,logfile,runNow,runCluster);
    pause(0.2);
end

%% Collect and Plot results.
% (Do this after the rest has completed)


Pars = zeros(length(geneNames),9);
Fano1 = zeros(1,length(geneNames));
Mean1 = zeros(1,length(geneNames));
Var1 = zeros(1,length(geneNames));
FanoEnd = zeros(1,length(geneNames));
MeanEnd = zeros(1,length(geneNames));
VarEnd = zeros(1,length(geneNames));
DeltaFano = zeros(1,length(geneNames));
DeltaMean = zeros(1,length(geneNames));
DeltaVar = zeros(1,length(geneNames));
for iGene = 1:length(geneNames)
    modelName = ['Model_',geneNames{iGene}];
    saveName = ['seqModels/',modelName];
    % load(saveName)
    eval(['Pars(iGene,:)=[',modelName,'.parameters{:,2}];'])
    Mean1(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.mean(1)']);
    Var1(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(1)']);
    Fano1(iGene) = Var1(iGene)/Mean1(iGene);
    CV1(iGene) = Var1(iGene)/Mean1(iGene)^2;
    
    MeanEnd(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.mean(end)']);
    VarEnd(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)']);
    FanoEnd(iGene) = VarEnd(iGene)/MeanEnd(iGene);
    CVEnd(iGene) = VarEnd(iGene)/MeanEnd(iGene)^2;
    
    DeltaFano(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)/',...
        'Model_',geneNames{iGene},'.dataSet.mean(end) -',... 
        'Model_',geneNames{iGene},'.dataSet.var(1)/',...
        'Model_',geneNames{iGene},'.dataSet.mean(1)']);
    DeltaCV(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)/',...
        'Model_',geneNames{iGene},'.dataSet.mean(end)^2 -',... 
        'Model_',geneNames{iGene},'.dataSet.var(1)/',...
        'Model_',geneNames{iGene},'.dataSet.mean(1)^2']);
    DeltaMean(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.mean(end)/ ',...
         'Model_',geneNames{iGene},'.dataSet.mean(1)']);
    DeltaVar(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)/ ',...
         'Model_',geneNames{iGene},'.dataSet.var(1)']);
    
end

% plots = {'Mean1','Var1','Fano1','MeanEnd','VarEnd','FanoEnd','DeltaMean','DeltaVar','DeltaFano'};
% plots = {'Mean1','Var1','CV1','MeanEnd','VarEnd','CVEnd','DeltaMean','DeltaVar','DeltaCV'};
% 
% figure(1);clf;
% for i=1:9
%     subplot(3,3,i)
%     if i<=3
%         scatter(Pars(:,3),Pars(:,5),100,log10(eval(plots{i})),'filled');
%         xlabel('K_{ON}'); ylabel('K_{OFF}');
%     elseif i<=6
%         scatter(Pars(:,3)+Pars(:,4),Pars(:,5)./(1+Pars(:,6)),100,log10(eval(plots{i})),'filled');
%         xlabel('K_{ON}'); ylabel('K_{OFF}');
%     else
%         scatter(Pars(:,4)./Pars(:,3),Pars(:,6),100,log10(eval(plots{i})),'filled');
%         xlabel('K_{ON}^{(1)}/K_{ON}^{(0)}');
%         ylabel('\alpha_{OFF}');
%     end
% end  
% subplot(3,3,1)
% title('Expression Mean')
% 
% subplot(3,3,2)
% title('Expression Variance')
% 
% subplot(3,3,3)
% title('Fano Factor')
% 
% for i=1:9
%     subplot(3,3,i)
%     set(gca,'XScale','log','YScale','log','FontSize',16)
%     if i<=6
%         set(gca,'xlim',[1e-5,1e2],'ylim',[1e-5,1e5])
%     elseif i>=7
%         xlim = get(gca,'XLim');ylim = get(gca,'YLim'); hold on
%         plot([1,1],ylim,'k--',xlim,[1,1],'k--','LineWidth',2);
%     end
%     colorbar
% end

plots = {'CV1','CVEnd','DeltaCV'};

figure(1);clf;
for i=1:3
    subplot(3,1,i)
    if i==1
        scatter(Pars(:,3),Pars(:,5),100,log10(eval(plots{i})),'filled');
        xlabel('K_{ON}'); ylabel('K_{OFF}');
    elseif i==2
        scatter(Pars(:,3)+Pars(:,4),Pars(:,5)./(1+Pars(:,6)),100,log10(eval(plots{i})),'filled');
        xlabel('K_{ON}'); ylabel('K_{OFF}');
    else
        scatter(Pars(:,4)./Pars(:,3),1./(1+Pars(:,6)),100,log10(eval(plots{i})),'filled');
        xlabel('K_{ON}^{(1)}/K_{ON}^{(0)}');
        ylabel('K_{OFF}^{(1)}/K_{OFF}^{(0)}');
        % ylabel('\alpha_{OFF}');
    end
end  
subplot(3,1,1)
title('Initial Time')

subplot(3,1,2)
title('Final Time')

subplot(3,1,3)
title('Change')

for i=1:3
    subplot(3,1,i)
    set(gca,'XScale','log','YScale','log','FontSize',16)
    if i<=2
        set(gca,'xlim',[1e-5,1e2],'ylim',[1e-5,1e5])
    else
        xlim = get(gca,'XLim');ylim = get(gca,'YLim'); hold on
        plot([2,2],ylim,'k--',xlim,[1/2,1/2],'k--','LineWidth',2);
    end
    C = colorbar;
    C.Label.String = 'CV^2';

end

%% plot input signal.
r1 = Model_ZNF689.parameters{1,2};
r2 = Model_ZNF689.parameters{2,2};
t = linspace(0,18,100);
Input = exp(-r1*t).*(1-exp(-r2*t));
figure(2)
plot(t,Input)
