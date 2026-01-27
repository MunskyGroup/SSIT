function makeMarginalsFsp(app)
% This function creates graphs of the different marginals at the desired
% time point in the FSP 3 Species GUI. This will use the time selected from
% the time slider as the time point from each of the different combinations
% based on the species selected.

T_array = eval(app.FspPrintTimesField.Value);
[~,j] = min(abs(T_array-app.FspTimeSlider.Value));

% Compute the marginal distributions
Nd = app.SSITModel.Solutions.fsp{j}.p.dim;
if Nd==1
    mdist{1} = double(app.SSITModel.Solutions.fsp{j}.p.data);
else
    for iSpecies=1:Nd
        INDS = setdiff([1:Nd],iSpecies);
        if ~isempty(INDS)
            mdist{iSpecies} = double(app.SSITModel.Solutions.fsp{j}.p.sumOver(INDS).data);
        else
            mdist{iSpecies} = double(app.SSITModel.Solutions.fsp{j}.p.data);
        end
    end
end

%% Check What species to plot
species2Plot=[];
legends={};
nSpecies = length(app.SpeciestoShowListBoxMargFSP.Items);
for iSpecies = 1:nSpecies
    if max(strcmpi(app.SpeciestoShowListBoxMargFSP.Value,app.SpeciestoShowListBoxMargFSP.Items{iSpecies}))
        species2Plot = [species2Plot iSpecies];
        legends=[legends char(app.SpeciestoShowListBoxMargFSP.Items{iSpecies})];
    end
end

speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);

% If needed compute distorted distributions.
if isfield(app.SSITModel.pdoOptions,'PDO')&&max(species2Plot)>length(speciesStochastic)
    maxNum = app.SSITModel.Solutions.fsp{j}.p.data.size;
    kSp = 0;
    for iS = 1:length(speciesStochastic)
        if max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic(iS)))
            maxNum(iS) = 0;
            curNum(iS) = 0;
        else
            kSp = kSp+1;
            curNum(iS) = size(app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp},2);
        end
    end
    if max(maxNum-curNum)>0
        [~,app.SSITModel] = app.SSITModel.generatePDO([],[],[],[],maxNum);
    end

    kSp = 0;
    for iSp = 1:length(speciesStochastic)
        if ~max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic{iSp}))
            INDS = setdiff([1:Nd],iSp);
            if ~isempty(INDS)
                px = double(app.SSITModel.Solutions.fsp{j}.p.sumOver(INDS).data);
            else
                px = double(app.SSITModel.Solutions.fsp{j}.p.data);
            end
            kSp = kSp+1;
            mdist{end+1} = app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp}*px;
        end
    end
end

figure()
hold('off')
cols = {'b', 'r', 'g', 'm', 'c', 'b--', 'r--', 'g--', 'm--', 'c--'};
legs = {};
for iSpecies = species2Plot
    stairs([0:length(mdist{iSpecies})], [mdist{iSpecies};0],cols{iSpecies},'LineWidth',2);  hold('on');
    legs{end+1} = app.SpeciestoShowListBoxMargFSP.Items{iSpecies};
end
title(sprintf('Marginals at time: t = %1.2f',T_array(j)));
xlabel('Species Count')
ylabel('Probability')
legend(legs);

set(gca,'fontsize',20);
hold('off');

end