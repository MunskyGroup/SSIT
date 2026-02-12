%Multi-Model Example.
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


%% Call Pipeline to Fit Model
% See ClusteScriptOnly.m

%% Collect and Plot results
% (Do this after the rest has completed)
% clear
% addpath(genpath('../src'));

DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
TAB = readtable(DataFileName);
geneNames = fields(TAB);
geneNames = geneNames(2:end-4); % The Genes are in Columns 2 -> N-4

Pars = zeros(length(geneNames),9);
Fano1 = zeros(1,length(geneNames));
Mean1 = zeros(1,length(geneNames));
Var1 = zeros(1,length(geneNames));
CV1 = zeros(1,length(geneNames));
FanoEnd = zeros(1,length(geneNames));
CVEnd = zeros(1,length(geneNames));
MeanEnd = zeros(1,length(geneNames));
VarEnd = zeros(1,length(geneNames));
DeltaFano = zeros(1,length(geneNames));
DeltaMean = zeros(1,length(geneNames));
MaxDeltaMean = zeros(1,length(geneNames));
DeltaVar = zeros(1,length(geneNames));
DeltaCV = zeros(1,length(geneNames));
fileExists = zeros(1,length(geneNames));
for iGene = 1:length(geneNames)
    modelName = ['Model_',geneNames{iGene}];
    saveName = ['seqModels/',modelName,'.mat'];
    if exist(saveName,'file')
        fileExists(iGene) = true;
        load(saveName)
        eval(['Pars(iGene,:)=[',modelName,'.parameters{:,2}];'])
        Mean1(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.mean(1)'])+1e-6;
        Var1(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(1)']);
        Fano1(iGene) = Var1(iGene)/Mean1(iGene);
        CV1(iGene) = Var1(iGene)/Mean1(iGene)^2;

        MeanEnd(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.mean(end)'])+1e-6;
        VarEnd(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)']);
        FanoEnd(iGene) = VarEnd(iGene)/MeanEnd(iGene);
        CVEnd(iGene) = VarEnd(iGene)/MeanEnd(iGene)^2;

        DeltaFano(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)/',...
            'Model_',geneNames{iGene},'.dataSet.mean(end) -',...
            'Model_',geneNames{iGene},'.dataSet.var(1)/',...
            'Model_',geneNames{iGene},'.dataSet.mean(1)']);
        MaxDeltaFano(iGene) = eval(['max(Model_',geneNames{iGene},'.dataSet.var./',...
            'Model_',geneNames{iGene},'.dataSet.mean) -',...
            'Model_',geneNames{iGene},'.dataSet.var(1)/',...
            'Model_',geneNames{iGene},'.dataSet.mean(1)']);
        DeltaCV(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)/',...
            'Model_',geneNames{iGene},'.dataSet.mean(end)^2 -',...
            'Model_',geneNames{iGene},'.dataSet.var(1)/',...
            'Model_',geneNames{iGene},'.dataSet.mean(1)^2']);
        DeltaMean(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.mean(end)/ ',...
            'Model_',geneNames{iGene},'.dataSet.mean(1)']);
        MaxDeltaMean(iGene) = eval(['max(Model_',geneNames{iGene},'.dataSet.mean)/ ',...
            'Model_',geneNames{iGene},'.dataSet.mean(1)']);
        DeltaVar(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)/ ',...
            'Model_',geneNames{iGene},'.dataSet.var(1)']);
        DeltaVar(iGene) = eval(['Model_',geneNames{iGene},'.dataSet.var(end)/ ',...
            'Model_',geneNames{iGene},'.dataSet.var(1)']);
    else
        fileExists(iGene) = false;
    end

end

%% Make figures
Xvals = (Pars(:,3) + Pars(:,4)*0.956)./Pars(:,3);  % activation ratio KON;
Yvals = (1 + Pars(:,6)*0.956); % activation ratio KOFF

Xvals1 = Pars(:,3);
Xvals2 = (Pars(:,3) + Pars(:,4)*0.956);

Yvals1 = Pars(:,5);
Yvals2 = Pars(:,5)./(1 + Pars(:,6)*0.956);

figure(1); clf;
J = 0*Xvals1;
for i = 1:length(Xvals1)
    if (Xvals2(i)/Xvals1(i)) >= 3*(Yvals1(i)/Yvals2(i))
        subplot(2,2,3);
        plot([Xvals1(i),Xvals2(i)],[Yvals1(i),Yvals2(i)],'-s'); hold on
        J(i) = 1;
    elseif (Xvals2(i)/Xvals1(i)) <= 1/3*(Yvals1(i)/Yvals2(i))
        subplot(2,2,2);
        plot([Xvals1(i),Xvals2(i)],[Yvals1(i),Yvals2(i)],'-o'); hold on
        J(i) = 2;
    else 
        subplot(2,2,4);
        plot([Xvals1(i),Xvals2(i)],[Yvals1(i),Yvals2(i)],'-^'); hold on
        J(i) = 3;
    end
    set(gca,'xscale','log','yscale','log','FontSize',16,...
        'ylim',[1e-4,1e4],'xlim',[1e-6,1e2],...
        'XTick',10.^[-6:4:2],'YTick',10.^[-4:4:4])
    xlabel('K_{ON}')
    ylabel('K_{OFF}')

end

subplot(2,2,1); hold off
scatter(Xvals(J==1),Yvals(J==1),100,log2(MaxDeltaMean(J==1))+1e-6,'s','filled'); hold on;
scatter(Xvals(J==2),Yvals(J==2),100,log2(MaxDeltaMean(J==2))+1e-6,'o','filled'); hold on;
scatter(Xvals(J==3),Yvals(J==3),100,log2(MaxDeltaMean(J==3))+1e-6,'^','filled'); hold on;
set(gca,'xscale','log','yscale','log','FontSize',16,...
    'xtick',10.^[0:2:6],'xlim',[1e0,1e6])
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
plot(xlim,[2,2],'k--')
plot([2,2],ylim,'k--')
C = colorbar;
C.Limits = [0,9];
C.Ticks = 0:2:8;
C.TickLabels = {'2^0', '2^2', '2^4', '2^6', '2^8'};
C.Label.String= 'Mean Fold Change';
xlabel('max(K_{ON}) / min(K_{ON})')
ylabel('max(K_{OFF}) / min(K_{OFF})')


%% Ignore below
% plots = {'Mean1','Var1','Fano1','MeanEnd','VarEnd','FanoEnd','DeltaMean','DeltaVar','DeltaFano'};
% % plots = {'Mean1','Var1','CV1','MeanEnd','VarEnd','CVEnd','DeltaMean','DeltaVar','DeltaCV'};
% %
% figure(1);clf;
% for i=1:9
%     subplot(3,3,i)
%     if i<=3
%         scatter(Pars(:,3)/Pars(9),Pars(:,5)/Pars(9),100,log10(eval(plots{i})+1e-6),'filled');
%         xlabel('K_{ON}'); ylabel('K_{OFF}');
%     elseif i<=6
%         scatter((Pars(:,3)+Pars(:,4))/Pars(9),(Pars(:,5)./(1+Pars(:,6)))/Pars(9),100,log10(eval(plots{i})),'filled');
%         xlabel('K_{ON}'); ylabel('K_{OFF}');
%     else
%         scatter(Pars(:,4)./Pars(:,3),Pars(:,6),100,log10(eval(plots{i})+1e-6),'filled');
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
%         set(gca,'xlim',[1e-5,1e3],'xtick',10.^[-5:2:3],'ylim',[1e-5,1e6],'ytick',10.^[-5:3:6])
%     elseif i>=7
%         xlim = get(gca,'XLim');ylim = get(gca,'YLim'); hold on
%         plot([1,1],ylim,'k--',xlim,[1,1],'k--','LineWidth',2);
%     end
%     colorbar
% end
% %%
% plots = {'CV1','CVEnd','DeltaCV'};
% 
% figure(1);clf;
% for i=1:3
%     subplot(3,1,i)
%     if i==1
%         scatter(Pars(:,3),Pars(:,5),100,log10(eval(plots{i}))','filled');
%         xlabel('K_{ON}'); ylabel('K_{OFF}');
%     elseif i==2
%         scatter(Pars(:,3)+Pars(:,4),Pars(:,5)./(1+Pars(:,6)),100,log10(eval(plots{i}))','filled');
%         xlabel('K_{ON}'); ylabel('K_{OFF}');
%     else
%         scatter((Pars(:,3)+Pars(:,4))./Pars(:,3),1./(1+Pars(:,6)),100,log10(eval(plots{i}))','filled');
%         xlabel('K_{ON}^{(1)}/K_{ON}^{(0)}');
%         ylabel('K_{OFF}^{(1)}/K_{OFF}^{(0)}');
%         % ylabel('\alpha_{OFF}');
%     end
% end
% subplot(3,1,1)
% title('Initial Time')
% 
% subplot(3,1,2)
% title('Final Time')
% 
% subplot(3,1,3)
% title('Change')
% 
% for i=1:3
%     subplot(3,1,i)
%     set(gca,'XScale','log','YScale','log','FontSize',16)
%     if i<=2
%         set(gca,'xlim',[1e-5,1e2],'ylim',[1e-5,1e5])
%     else
%         xlim = get(gca,'XLim');ylim = get(gca,'YLim'); hold on
%         plot([2,2],ylim,'k--',xlim,[1/2,1/2],'k--','LineWidth',2);
%     end
%     C = colorbar;
%     C.Label.String = 'CV^2';
% 
% end
% 
% %% plot input signal.
% r1 = Model_ZNF689.parameters{1,2};
% r2 = Model_ZNF689.parameters{2,2};
% t = linspace(0,18,100);
% Input = exp(-r1*t).*(1-exp(-r2*t));
% figure(2)
% plot(t,Input)
