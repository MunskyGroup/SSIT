function makePlotsDUSP1(ModelDusp1Fit,ModelGRDusp,DUSP1pars,Dusp1FitCases,showCases)
if showCases(1)
    fignums = [211,221,201,231];
    % combinedDusp1Model = combinedDusp1Model.updateModels(DUSP1pars,false,fignums);
    for i=1:size(Dusp1FitCases,1)
        %  Update parameters in original models.
        ModelDusp1Fit{i}.parameters(ModelDusp1Fit{i}.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        ModelDusp1Fit{i}.tSpan = sort(unique([ModelDusp1Fit{i}.tSpan,linspace(0,180,30)]));

        ModelDusp1Fit{i}.makeFitPlot([],5,fignums)
        figure(Dusp1FitCases{i,3});
        set(gca,'ylim',[0,150])
        title(Dusp1FitCases{i,4})
        ylabel('Nuclear DUSP1 mRNA')
        xlabel('Time (min)')
    end
end

%%  PREDICT DUSP1 Distributions under other Dex concentrations.
if showCases(2)
    PredictionCases = {'10','10',301,'DUSP1 Prediction (10nM Dex)';...
        '1','1',302,'DUSP1 Prediction (1nM Dex)';...
        '0p3','0.3',303,'DUSP1 Prediction (0.3nM Dex)'};

    fignums = [311,321,301,331;...
        312,322,302,332;...
        313,323,303,333];
    ModelPred = cell(size(Dusp1FitCases,1),1);
    for i=1:size(PredictionCases,1)
        ModelPred{i} = ModelGRDusp.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
            {'rna','totalNucRNA'},...
            {'Dex_Conc',PredictionCases{i,2}});

        ModelPred{i}.tSpan = sort(unique([ModelPred{i}.tSpan,linspace(0,180,30)]));

        if str2num(PredictionCases{i,2})~=100
            % ModelPred{i}.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor = ...
            NT = size(ModelPred{i}.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor,1)+1;
            NS = size(ModelPred{i}.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor,2);
            TMP(1,:) = double(ModelGRDusp.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(1,:));
            TMP(2:NT,1:NS) = double(ModelPred{i}.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor); 
            ModelPred{i}.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor = sptensor(TMP);
            ModelPred{i}.dataSet.times = unique([0,ModelPred{i}.dataSet.times]);
        end

        % Set model parameters to those supplied
        ModelPred{i}.parameters(ModelPred{i}.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        % Change Dex concentraion in the model.
        ModelPred{i}.parameters(13,1:2) = {'Dex0',str2num(PredictionCases{i,2})};

        ModelPred{i}.makeFitPlot([],5,fignums(i,:))

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

    ModelPredDexTtrSoln = cell(size(DexConc,1),1);
    for i=1:length(DexConc)
        ModelPredDexTtr = ModelGRDusp;
        ModelPredDexTtr.dataSet = [];
        ModelPredDexTtr = ModelPredDexTtr.loadData('EricData/pdoCalibrationData_EricIntensity_ConcSweeps.csv',...
            {'rna','totalNucRNA'},...
            {'Dex_Conc',DecConcStr{i}});

        % Set model parameters to those supplied
        ModelPredDexTtr.parameters(ModelPredDexTtr.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        
        % Change Dex concentration in the model.
        ModelPredDexTtr.parameters(13,1:2) = {'Dex0',str2num(DecConcStr{i})};

        % ModelPredDexTtr = ModelPredDexTtr.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_TtrPred']);
        ModelPredDexTtrSoln = ModelPredDexTtr.solve;

        DataHist = double(ModelPredDexTtr.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor);
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
    fignums = [511,521,501,531;...
        512,522,502,532;...
        513,523,503,533;...
        514,524,504,534;...
        515,525,505,535];

    % List of tryptolide experiments
    PredictionCases = {'0',501,'DUSP1 Prediction (t_{TPL} = 0 min)';...
        '20',502,'DUSP1 Prediction (t_{TPL} = 20 min)';...
        '75',503,'DUSP1 Prediction (t_{TPL} = 75 min)';...
        '150',504,'DUSP1 Prediction (t_{TPL} = 150 min)';...
        '180',505,'DUSP1 Prediction (t_{TPL} = 180 min)'};

    ModelPredDexTpl = cell(size(PredictionCases,1),1);
    ModelPredDexTplSoln = cell(size(PredictionCases,1),1);
    for i=1:size(PredictionCases,1)
        ModelPredDexTpl{i}.dataSet = [];
        ModelPredDexTpl{i} = ModelGRDusp.loadData('EricData/DUSP1_predict_data_TPL_SSIT.csv',...
            {'rna','totalNucRNA'},...
            {['tryptCond',num2str(i)],num2str(i)});

        % Set model parameters to those supplied
        ModelPredDexTpl{i}.parameters(ModelPredDexTpl{i}.fittingOptions.modelVarsToFit,2) = num2cell(DUSP1pars);
        
        % set the Dex concentration.
        ModelPredDexTpl{i}.parameters(13,:) = {'Dex0',100};

        ModelPredDexTpl{i}.tSpan = sort(unique([ModelPredDexTpl{i}.tSpan,linspace(0,250,30)]));

        ModelPredDexTpl{i}.propensityFunctions = {'kon*offGene*nucGR';'koff*onGene';
            '(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';'kg1';'gg1*cytGR';'gg2*nucGR';...
            'kr*onGene*Itrpt';'gr*rna'};

        ModelPredDexTpl{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)';...
            'Itrpt',['(t<=',PredictionCases{i},')']};

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