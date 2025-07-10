function [KS] = makeCytDistPlots(ssaSoln,extendedMod,fignum,timeIndsDat,speciesIndMod,speciesIndDat,binEdges,splitReps)
arguments
    ssaSoln
    extendedMod
    fignum = 1;
    timeIndsDat = [];
    speciesIndMod = 1;
    speciesIndDat = 1;
    binEdges = [0:10:300];
    splitReps = false
end
figure(fignum); clf;
if isempty(timeIndsDat)
    timeIndsDat = [1:length(extendedMod.dataSet.times)];
end

times2plot = extendedMod.dataSet.times(timeIndsDat);

Nrows = ceil(sqrt(length(timeIndsDat)));
Ncols = ceil(length(timeIndsDat)/Nrows);
KS = zeros(1,length(timeIndsDat));
for i = 1:length(timeIndsDat)
    subplot(Ncols,Nrows,i); hold off
    time = times2plot(i);
    
    ctime = 13;
    creplica = 15;
    crnadat = [3,4];
    if splitReps
        dMat = extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,[crnadat,creplica]);
        repNames = unique(dMat(:,3));
        for j = 1:length(repNames)
            dMatB = cell2mat(dMat(strcmp(dMat(:,3),repNames{j}),1:2));
            h = histogram(dMatB(:,speciesIndDat),binEdges,'Normalization','pdf','DisplayStyle','stairs','LineWidth',2);
            hold on
        end
    else
        dMatB = cell2mat(extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,crnadat));
        histogram(dMatB(:,speciesIndDat),binEdges,'Normalization','pdf','DisplayStyle','stairs','LineWidth',3,'EdgeColor',[.40,.40,.80]);
        hold on
    end
    dMatB = cell2mat(extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,crnadat));

    % Add SSA to histogram plot
    % Find time in SSA data
    [~,jSp] = min(abs(time-ssaSoln.T_array));
    M = squeeze(ssaSoln.trajs(speciesIndMod,jSp,:));
    histogram(M,binEdges,'Normalization','pdf','DisplayStyle','stairs','LineWidth',3,'EdgeColor',[0,0,0]);

    [~,~,KS(i)] = kstest2(dMatB(:,speciesIndDat),M');
    % % Add data to histogram plot
    % dMat = double(extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(timeIndsDat(i),:,:));
    % N = sum(dMat,"all");
    % if speciesIndDat==1
    %     PD = [0;sum(dMat,2)/N];
    % else
    %     PD = [0;sum(dMat,1)'/N];
    % end
    % binEdges = round(H.BinEdges);
    % PD(binEdges(end)+1) = 0;
    % 
    % nBins = length(binEdges)-1;
    % PDbinned = zeros(nBins+1,1);
    % binwidth = binEdges(2)-binEdges(1);
    % for j = 1:nBins
    %     PDbinned(j) = sum(PD(binEdges(j)+1:binEdges(j+1)));
    % end
    % 
    % PDbinned(end) = 1-sum(PDbinned);
    % stairs(binEdges,PDbinned/binwidth,'linewidth',2)
    title(['t = ',num2str(time),'min']);

    set(gca,'FontSize',15,'ylim',[0,0.03])

    
end
end
