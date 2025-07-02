function [sumKS] = compareDistPlots(ssaSoln,extendedMod,timeIndsDat,speciesIndMod,speciesIndDat)
arguments
    ssaSoln
    extendedMod
    timeIndsDat = [];
    speciesIndMod = 1;
    speciesIndDat = 1;
end
if isempty(timeIndsDat)
    timeIndsDat = [1:length(extendedMod.dataSet.times)];
end

times2plot = extendedMod.dataSet.times(timeIndsDat);

KS = zeros(1,length(timeIndsDat));
for i = 1:length(timeIndsDat)
    time = times2plot(i);
    [~,jSp] = min(abs(time-ssaSoln.T_array));
    
    ctime = 13;
    crnadat = [3,4];
 
    dMatB = cell2mat(extendedMod.dataSet.DATA([extendedMod.dataSet.DATA{:,ctime}]==time,crnadat));

    M = squeeze(ssaSoln.trajs(speciesIndMod,jSp,:));
    [~,~,KSi] = kstest2(dMatB(:,speciesIndDat),M');
    KS(i) = KSi*length(dMatB);
    
end
sumKS = sum(KS);
end
