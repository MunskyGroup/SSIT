function cleanUpFile(filenameTMP)
load([filenameTMP,'.mat']);
b=whos;
varNames = {b.name};
varNames = varNames(~contains(varNames,'filename'));
save([filenameTMP,'_bak.mat'],varNames{:});
save([filenameTMP,'.mat'],varNames{:});
end