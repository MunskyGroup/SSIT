function makeMeanAndVarianceFsp(app)
% This function creates graphs of the different marginals at the desired
% time point in the FSP 3 Species GUI. This will use the time selected from
% the time slider as the time point from each of the different combinations
% based on the species selected.

T_array = eval(app.FspPrintTimesField.Value);

%% Check What species to plot
species2Plot=[];
legends={};
nSpecies = length(app.SpeciestoShowListBoxMeans.Items);
for iSpecies = 1:nSpecies
    if max(strcmpi(app.SpeciestoShowListBoxMeans.Value,app.SpeciestoShowListBoxMeans.Items{iSpecies}))
        species2Plot = [species2Plot iSpecies];
        legends=[legends char(app.SpeciestoShowListBoxMeans.Items{iSpecies})];
    end
end
speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);

Nd = app.SSITModel.Solutions.fsp{1}.p.dim;
mdist = cell(length(T_array),nSpecies);
for it = 1:length(T_array)
    if Nd==1
        mdist{it,1} = double(app.SSITModel.Solutions.fsp{it}.p.data);
    else
        for iSpecies=1:Nd
            INDS = setdiff([1:Nd],iSpecies);
            if ~isempty(INDS)
                mdist{it,iSpecies} = double(app.SSITModel.Solutions.fsp{it}.p.sumOver(INDS).data);
            else
                mdist{it,iSpecies} = double(app.SSITModel.Solutions.fsp{it}.p.data);
            end
        end
    end

    % If needed compute distorted distributions.
    if isfield(app.SSITModel.pdoOptions,'PDO')&&max(species2Plot)>length(speciesStochastic)
        maxNum = app.SSITModel.Solutions.fsp{end}.p.data.size;
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
            if isempty(app.SSITModel.pdoOptions.unobservedSpecies)||~max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic{iSp}))
                INDS = setdiff([1:Nd],iSp);
                if ~isempty(INDS)
                    px = double(app.SSITModel.Solutions.fsp{it}.p.sumOver(INDS).data);
                else
                    px = double(app.SSITModel.Solutions.fsp{it}.p.data);
                end
                kSp = kSp+1;
                mdist{it,Nd+kSp} = app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp}(:,1:length(px))*px;
            end
        end
    end

end

for it = 1:length(T_array)
    for j=1:max(species2Plot)
        mns(it,j) = [0:length(mdist{it,j})-1]*mdist{it,j};
        mns2(it,j) = [0:length(mdist{it,j})-1].^2*mdist{it,j};
    end
end
vars = mns2-mns.^2;


%% Get ODE solutions if needed.
if app.FspMeanVarShowOdeCheckBox.Value == 1
    if isempty(app.SSITModel.Solutions)||~isfield(app.SSITModel.Solutions,'ode')
        currentSolnScheme = app.SSITModel.solutionScheme;
        app.SSITModel.solutionScheme = 'ode';
        [~,~,app.SSITModel] = app.SSITModel.solve;
        app.SSITModel.solutionScheme = currentSolnScheme;
    end
    [~,iA] = intersect(app.SSITModel.species,speciesStochastic);
    odeSolns = app.SSITModel.Solutions.ode(:,iA);
end

figure();
cols = {'b','r','g','m','c','y','k'};
cols2 = {'--m','--c','--k','b--','r--','g--','m--','c--','y--','k--'};
LG = {};
for iplt=species2Plot
    if app.FspMeanVarShowVarianceCheckBox.Value == 1
        shadedErrorBar(T_array,mns(:,iplt),sqrt(vars(:,iplt)),'lineprops',cols{iplt});hold on;
        LG{end+1} = [app.SpeciestoShowListBoxMeans.Items{iplt},' \pm std'];
    else
        plot(T_array,mns(:,iplt),cols{iplt});hold on;
        LG{end+1} = app.SpeciestoShowListBoxMeans.Items{iplt};
    end
    if app.FspMeanVarShowOdeCheckBox.Value == 1&&iplt<=Nd
        % if isempty(app.FspTabOutputs.odeSolutions)
        plot(app.SSITModel.tSpan,odeSolns(:,iplt),cols2{iplt},'linewidth',2); hold on;
        plot(T_array,mns(:,iplt),cols{iplt});
        LG{end+1} = [append(app.SSITModel.species{iplt},' ode')];
        LG{end+1} = app.SSITModel.species{iplt};
    end

end
set(gca,'fontsize',20)
title('Averages Response vs Time');
xlabel('Time');
ylabel('Average Response');
legend(LG)

end
