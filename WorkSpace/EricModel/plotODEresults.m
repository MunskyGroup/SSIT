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
    plot(modeWithGRData.dataSet.times,modeWithGRData.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
    legend(extendedMod.species(3:4))
    set(gca,'xlim',[-10,200],'ylim',[0,12],'fontsize',16)
    ylabel('GR Concentrations (UA)')
    legend({'Cyt-Model','Nuc-Model','Cyt-Data','Nuc-Data'})
    title('GR')

    % Plot DUSP1 levels vs. Time
    subplot(2,1,2)
    plot(extendedMod.tSpan,soln.ode(:,5:6),'--','LineWidth',2);hold on
    plot(extendedMod.dataSet.times,extendedMod.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
    legend(extendedMod.species(3:4))
    set(gca,'xlim',[-10,200],'ylim',[0,160],'fontsize',16)
    ylabel('GR Concentrations (UA)')
    legend({'Nuc-Model','Cyt-Model','Nuc-Data','Cyt-Data'})
    title('DUSP1')
end