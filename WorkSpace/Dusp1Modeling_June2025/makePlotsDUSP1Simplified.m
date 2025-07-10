function makePlotsDUSP1Simplified(ModelDusp1Fit,DUSP1pars,showCases,kTPL,modelLibrary)
arguments
    ModelDusp1Fit
    DUSP1pars
    showCases
    kTPL = 0.693;
    modelLibrary = 'savedParameters/GRDusp1ModelLibrary';
end

if showCases(1)
    fignums = [211,221,201,231];
    % combinedDusp1Model = combinedDusp1Model.updateModels(DUSP1pars,false,fignums);
    %  Update parameters in original models.
    ModelDusp1Fit.parameters(ModelDusp1Fit.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
    ModelDusp1Fit.tSpan = sort(unique([ModelDusp1Fit.tSpan,linspace(0,180,30)]));

    ModelDusp1Fit.makeFitPlot([],5,fignums)
    figure(201);
    set(gca,'ylim',[0,150])
    title('DUSP1 Fit (100nM Dex)')
    ylabel('Nuclear DUSP1 mRNA')
    xlabel('Time (min)')
end

%%  PREDICT DUSP1 Distributions under other Dex concentrations.
if showCases(2)
    PredictionCases = {'10','10',301,'DUSP1 Prediction (10nM Dex)','ModelDUSP1_10nM';...
        '1','1',302,'DUSP1 Prediction (1nM Dex)','ModelDUSP1_1nM';...
        '0p3','0.3',303,'DUSP1 Prediction (0.3nM Dex)','ModelDUSP1_0p3nM'};

    fignums = [311,321,301,331;...
        312,322,302,332;...
        313,323,303,333];
    
    for i=1:3
        Model = dusp1ModelLibrary(PredictionCases{i,5},false,modelLibrary);
        Model.tSpan = sort(unique([Model.tSpan,linspace(0,180,30)]));

        % Set model parameters to those supplied
        Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        
        % Make plots
        Model.makeFitPlot([],5,fignums(i,:))

        figure(fignums(i,3));
        set(gca,'ylim',[0,150])
        title(PredictionCases{i,4})
        ylabel('DUSP1 mRNA')
        xlabel('Time (min)')
    end
end

%%  Predict DUSP1 Distributions at 75 min under Dex Titration
if showCases(3)
    DexConc = 10.^[-3,-2,-1,0,1,2,3,4];
    DecConcStr = {'0.001','0.01','0.1','1','10','100','1000','10000'};
    ModelPredDexTtr = dusp1ModelLibrary('ModelPredDexTtr',false,modelLibrary);

    for i=1:length(DexConc)
        Model = ModelPredDexTtr{i};
        
        % Set model parameters to those supplied
        Model.parameters(Model.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        
        ModelPredDexTtrSoln = Model.solve;

        DataHist = double(Model.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor);
        PModel = double(ModelPredDexTtrSoln.fsp{2}.p.sumOver([1,2]).data);

        MeanData(i) = DataHist*[0:length(DataHist)-1]'/sum(DataHist);
        MeanModel(i) = PModel'*[0:length(PModel)-1]'/sum(PModel);

        Mean2Data = DataHist*([0:length(DataHist)-1]').^2/sum(DataHist);
        Mean2Model = PModel'*([0:length(PModel)-1]').^2/sum(PModel);

        SigDat(i) = sqrt(Mean2Data-MeanData(i)^2);
        SigMod(i) = sqrt(Mean2Model-MeanModel(i)^2);

    end
    figure(401); clf; hold on
    x = [DexConc,DexConc(end:-1:1)];
    y = [MeanModel-SigMod,MeanModel(end:-1:1)+SigMod(end:-1:1)];
    fill(x,y,[0.9,1,0.9])
    plot(DexConc,MeanModel,'k','LineWidth',3); hold on;

    errorbar(DexConc,MeanData,SigDat,'LineWidth',3,'LineStyle','none'); hold on;
    plot(DexConc,MeanData,'s','MarkerSize',18,'MarkerFaceColor','b');

    set(gca,'xscale','log','FontSize',16,'xlim',[9e-4,1.1e4],'ylim',[0,150])
    title('Prediction of DUSP1 Expression')
    xlabel('Dex Concentration at 75 min')
    ylabel('DUSP1 Expression (nM)')
end

%%  Predict DUSP1 Distributions After Tryptolide
if showCases(4)
    fignums = [2011,2021,2001,2031;...
        2012,2022,2002,2032;...
        2013,2023,2003,2033;...
        2014,2024,2004,2034;...
        2015,2025,2005,2035];

    % List of tryptolide experiments
    PredictionCases = {'0',2001,'DUSP1 Prediction (t_{TPL} = 0 min)';...
        '20',2002,'DUSP1 Prediction (t_{TPL} = 20 min)';...
        '75',2003,'DUSP1 Prediction (t_{TPL} = 75 min)';...
        '150',2004,'DUSP1 Prediction (t_{TPL} = 150 min)';...
        '180',2005,'DUSP1 Prediction (t_{TPL} = 180 min)'};

    ModelPredDexTpl = cell(size(PredictionCases,1),1);
    ModelPredDexTplSoln = cell(size(PredictionCases,1),1);
    for i=1:size(PredictionCases,1)
        ModelPredDexTpl{i} = ModelGRDusp;
        ModelPredDexTpl{i}.dataSet = [];
        ModelPredDexTpl{i} = ModelPredDexTpl{i}.loadData('EricData/TryptolideData.csv',...
            {'rna','RNA_DUSP1_nuc'},...
            {'Time_TPL',PredictionCases{i,1}});

        % Set model parameters to those supplied
        ModelPredDexTpl{i}.parameters(ModelPredDexTpl{i}.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        
        % set the Dex concentration.
        ModelPredDexTpl{i}.parameters(13,:) = {'Dex0',100};

        ModelPredDexTpl{i}.tSpan = sort(unique([ModelPredDexTpl{i}.tSpan,linspace(0,250,30)]));

        ModelPredDexTpl{i}.propensityFunctions{8} = 'kr*onGene*Itrpt';

        % ModelPredDexTpl{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)';...
        %     'Itrpt',['(t<=',PredictionCases{i},')']};
        ModelPredDexTpl{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)';...
            'Itrpt',['(t<=',PredictionCases{i},') + (t>',PredictionCases{i},')*exp(-',num2str(kTPL),'*(t-',PredictionCases{i},'))']};

        ModelPredDexTpl{i} = ModelPredDexTpl{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_TplPred'],false);

        ModelPredDexTpl{i}.makeFitPlot([],5,fignums(i,:))

        figure(fignums(i,3));
        set(gca,'ylim',[0,150])
        title(PredictionCases{i,3})
        ylabel('DUSP1 mRNA')
        xlabel('Time (min)')
    end
end

end