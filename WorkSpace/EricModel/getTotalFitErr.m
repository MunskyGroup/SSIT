function output = getTotalFitErr(Organization,pars,init)
arguments
    Organization
    pars = []
    init = false

end

if init
    for i=1:size(Organization,1)
        Organization{i,1}.fittingOptions.modelVarsToFit = Organization{i,2};
        Organization{i,1}.fittingOptions.logPrior = [];
    end
    output = Organization;
    return
else
    parfor i=1:size(Organization,1)
        if isempty(pars)
            parsLocal = [Organization{i,1}.parameters{Organization{i,2},2}];
        else
            parsLocal = pars(Organization{i,3});
        end
        newJ(i) =  Organization{i,5}*Organization{i,1}.(Organization{i,4})(parsLocal)
    end
    output = sum(newJ);
end
end
