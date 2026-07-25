function FIM = totalFim(fims,Nc,covPrior)
arguments
    fims
    Nc
    covPrior = [];
end
if isempty(covPrior)
    fimPrior = zeros(size(fims{1}));
else
    fimPrior = inv(covPrior);
end
Nt = size(fims,1);
Ns = size(fims,2);
FIM = cell(1,Ns);
for is = 1:Ns
    FIM{is} = fimPrior;
    for it = 1:Nt
        FIM{is} = FIM{is}+Nc(it)*fims{it,is};
    end
end

end
