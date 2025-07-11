function makeNucCytScatterPlots(ssaSoln,extendedMod,fignum,speciesIndMod,speciesIndDat,splitReps)
arguments
    ssaSoln
    extendedMod
    fignum = 1;
    speciesIndMod = [1,2];
    speciesIndDat = [1,2];
    splitReps = false;
end
figure(fignum); clf;
    timeIndsDat = [1:length(extendedMod.dataSet.times)];

times2plot = extendedMod.dataSet.times(timeIndsDat);

Nrows = ceil(sqrt(length(timeIndsDat)));
Ncols = ceil(length(timeIndsDat)/Nrows);
for i = 1:length(timeIndsDat)
    subplot(Ncols,Nrows,i)

    % Find time in SSA data
    time = times2plot(i);
    [~,jSp] = min(abs(time-ssaSoln.T_array));
    
    % Add SSA to histogram plot
    M = squeeze(ssaSoln.trajs(speciesIndMod,jSp,:));
    % scatter(M(1,:),M(2,:),'ro');
    ksdensity(M','PlotFcn','contour');
    hold on

    % Add data to histogram plot
    ctime = 13;
    creplica = 15;
    crnadat = [3,4];
    if splitReps
        dMat = extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,[crnadat,creplica]);
        repNames = unique(dMat(:,3));

        for j = 1:length(repNames)
            dMatB = cell2mat(dMat(strcmp(dMat(:,3),repNames{j}),1:2));
            scatter(dMatB(:,1),dMatB(:,2),80,'.');
            hold on
        end
    else
        dMatB = cell2mat(extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,crnadat));
        scatter(dMatB(:,1),dMatB(:,2),80,'k.');
        hold on
    end
    xyAll = cell2mat(extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,crnadat));

    % dMat = extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor;
    % 
    % xy = dMat.subs(dMat.subs(:,1)==timeIndsDat(i),1+speciesIndDat);   
    % v = dMat.values(dMat.subs(:,1)==timeIndsDat(i));
    % xyAll = zeros(sum(v),2);
    % for j = 1:length(v)
    %     xyAll(sum(v(1:j-1))+1:sum(v(1:j)),:) = repmat(xy(j,:),v(j),1);
    % end    
    % scatter(xyAll(:,1),xyAll(:,2),20,'k.');

    par0 = mean(M');
    cov12 = cov(M'); 
    ssit.parest.ellipse(par0,icdf('chi2',0.68,2)*cov12,'b','linewidth',3);  hold on;

    par0 = mean(xyAll);
    cov12 = cov(xyAll); 
    ssit.parest.ellipse(par0,icdf('chi2',0.68,2)*cov12,'c--','linewidth',3);  hold on;
    
    title(['t = ',num2str(time),'min']);
    xlabel('Nuc');ylabel('Cyt');
    set(gca,'xlim',[0,200],'ylim',[0,300])
    
end
end
