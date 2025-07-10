function [sumKS,KS] = compareDistPlots(ssaSoln,extendedMod,speciesIndMod,speciesIndDat,bins,crnadat,ctime,creplica)
arguments
    ssaSoln
    extendedMod
    speciesIndMod = 1;
    speciesIndDat = 1;
    bins = [0:5:300];
    crnadat = [3,4];
    ctime = [13];
    creplica = [15];
end
timeIndsDat = [1:length(extendedMod.dataSet.times)];

times2plot = extendedMod.dataSet.times(timeIndsDat);

KS = zeros(1,length(timeIndsDat));
for i = 1:length(timeIndsDat)
    time = times2plot(i);
    [~,jSp] = min(abs(time-ssaSoln.T_array));
        
    dMat = extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,[crnadat,creplica]);
    repNames = unique(dMat(:,3));
    for j = 1:length(repNames)
        dMatB = cell2mat(dMat(strcmp(dMat(:,3),repNames{j}),1:2));
        if j==1
            Hmax = cumsum(hist(dMatB(:,speciesIndDat),bins)); Hmax=Hmax/Hmax(end);
            Hmin = Hmax;
        else
            H= cumsum(hist(dMatB(:,speciesIndDat),bins)); H=H/H(end);
            Hmin = min(Hmin,H);
            Hmax = max(Hmax,H);
        end
        hold on
    end    

    % dMatB = cell2mat(extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,crnadat));

    M = squeeze(ssaSoln.trajs(speciesIndMod,jSp,:));

    H1 = cumsum(hist(M,bins)); H1=H1/H1(end);
    % H2 = cumsum(hist(dMatB(:,speciesIndDat),[0:5:300])); H2=H2/H2(end);

    AUC = sum(max(H1-Hmax,0)) + sum(max(Hmin-H1,0));

    % AUC = sum(abs(H1-H2))*length(dMatB);
    KS(i) = AUC;
    % [~,~,KSi] = kstest2(dMatB(:,speciesIndDat),M');
    % KS(i) = KSi*length(dMatB);
    
end
sumKS = sum(KS);
end
