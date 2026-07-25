function k = findBestMove(fims,Ncp,met,NcMax,statistic,covPrior,incrementAdd)
arguments
    fims
    Ncp
    met
    NcMax = []
    statistic='mean'
    covPrior=[]
    incrementAdd=1
end

if isempty(statistic)
    statistic='mean';
end
Nt = size(fims,1);
if isempty(NcMax)
    NcMax = inf*ones(1,Nt);
end

Ns = size(fims,2);
objFun = zeros(Nt,Ns);
FIM0 = SSIT.totalFim(fims,Ncp,covPrior);

for is = 1:Ns
    for it = 1:Nt
        if Ncp(it)+incrementAdd<=NcMax(it)
            % If one can do that experiment.
            FIM = FIM0{is}+incrementAdd*fims{it,is};
        else
            % If there are no more cells avalable for that time
            % point.
            FIM = FIM0{is};
        end
        objFun(it,is) = met(FIM);
    end
end

switch statistic
    case 'median'
        [~,k] = min(median(objFun,2));
    case 'mean'
        [~,k] = min(mean(objFun,2));
    otherwise
        [~,k] = min(median(objFun,2));
end

end