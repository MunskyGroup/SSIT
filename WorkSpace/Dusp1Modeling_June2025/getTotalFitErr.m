function output = getTotalFitErr(Organization,pars,init,extraObjs)
arguments
    Organization
    pars = []
    init = false
    extraObjs = [];
end

if init
    for i=1:size(Organization,1)
        Organization{i,1}.fittingOptions.modelVarsToFit = Organization{i,2};
        Organization{i,1}.fittingOptions.logPrior = [];
    end
    output = Organization;
    return
else
    nObj = size(Organization,1)+size(extraObjs,1);
    newJ = zeros(1,nObj);
    for i=1:nObj
        if i<=size(Organization,1)
            if isempty(pars)
                parsLocal = [Organization{i,1}.parameters{Organization{i,2},2}];
            else
                parsLocal = pars(Organization{i,3});
            end
            newJ(i) =  Organization{i,5}*Organization{i,1}.(Organization{i,4})(parsLocal);
        else
            parsLocal = pars(extraObjs{i-size(Organization,1),2});
            newJ(i) =  extraObjs{i-size(Organization,1),1}(parsLocal);
        end
    end
    output = sum(newJ);
end
end
