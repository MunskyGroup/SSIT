function plotODEresults(extendedMod,soln,modeWithGRData,fignum)
arguments
    extendedMod
    soln
    modeWithGRData
    fignum = 1;
end
    figure(fignum); clf;
    % Plot GR levels vs. Time
    subplot(2,1,1)
    plot(extendedMod.tSpan,soln.ode(:,3:4),'--','LineWidth',2);hold on
    % plot(modeWithGRData.dataSet.times,modeWithGRData.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
    reps = {'A','B','C'};
    for i = 1:length(modeWithGRData.dataSet.times)
        for j = 1:3
            dataSelected = modeWithGRData.dataSet.DATA([modeWithGRData.dataSet.DATA{:,5}]'==modeWithGRData.dataSet.times(i)& ...
                strcmp(modeWithGRData.dataSet.DATA(:,7),reps{j}),:);
            meanNuc(j,i) = mean([dataSelected{:,11}]);
            meanCyt(j,i) = mean([dataSelected{:,10}]);
        end
    end
    mnNucAll = mean(meanNuc);
    stdNuc = std(meanNuc);
    mnCytAll = mean(meanCyt);
    stdCyt = std(meanCyt);

    errorbar(modeWithGRData.dataSet.times,mnNucAll,stdNuc,'LineWidth',3)
    errorbar(modeWithGRData.dataSet.times,mnCytAll,stdCyt,'LineWidth',3)
    legend(extendedMod.species(3:4))
    set(gca,'xlim',[-10,200],'ylim',[0,20],'fontsize',16)
    ylabel('GR Concentrations (UA)')
    legend({'Cyt-Model','Nuc-Model','Cyt-Data','Nuc-Data'})
    title('GR')

    % Plot DUSP1 levels vs. Time
    for i = 1:length(extendedMod.dataSet.times)
        downSelectData = extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,13}]'==extendedMod.dataSet.times(i),:);
        reps = unique(downSelectData(:,15));

        for j = 1:length(reps)
            dataSelected = downSelectData(strcmp(downSelectData(:,15),reps{j}),:);
            meanNucDusp1(j,i) = mean([dataSelected{:,3}]);
            meanCytDusp1(j,i) = mean([dataSelected{:,4}]);
        end
    end
    meanNucDusp1(meanNucDusp1==0) = NaN;
    meanCytDusp1(meanCytDusp1==0) = NaN;
    
    mnNucAllDusp1 = mean(meanNucDusp1,"omitmissing");
    stdNucDusp1 = std(meanNucDusp1,"omitmissing");
    mnCytAllDusp1 = mean(meanCytDusp1,"omitmissing");
    stdCytDusp1 = std(meanCytDusp1,"omitmissing");

    subplot(2,1,2)
    plot(extendedMod.tSpan,soln.ode(:,5:6),'--','LineWidth',3);hold on
    % plot(extendedMod.dataSet.times,extendedMod.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
    % errorbar(extendedMod.dataSet.times,extendedMod.dataSet.mean,sqrt(extendedMod.dataSet.var),'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
    errorbar(extendedMod.dataSet.times,mnNucAllDusp1,stdNucDusp1,'s','MarkerSize',14,'MarkerFaceColor','k','LineWidth',2)
    errorbar(extendedMod.dataSet.times,mnCytAllDusp1,stdCytDusp1,'s','MarkerSize',14,'MarkerFaceColor','k','LineWidth',2)
    legend(extendedMod.species(3:4))
    set(gca,'xlim',[-10,200],'ylim',[0,200],'fontsize',16)
    ylabel('GR Concentrations (UA)')
    legend({'Nuc-Model','Cyt-Model','Nuc-Data','Cyt-Data'})
    title('DUSP1')
end