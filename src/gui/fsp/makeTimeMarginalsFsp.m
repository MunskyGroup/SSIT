function makeTimeMarginalsFsp(app)
% This function creates plots comparing the different marginals throughout
% time of a single species. ReactionsTabOutputs.inputs require vector format for the time
% points.

T_array = eval(app.FspPrintTimesField.Value);
T_array2 = eval(app.FspMarginalVecField.Value);

%% Check What species to plot
species2Plot=[];
legends={};
nSpecies = length(app.SpeciestoShowListBoxMargFSPvT.Items);
for iSpecies = 1:nSpecies
    if max(strcmpi(app.SpeciestoShowListBoxMargFSPvT.Value,app.SpeciestoShowListBoxMargFSPvT.Items{iSpecies}))
        species2Plot = [species2Plot iSpecies];
        legends=[legends char(app.SpeciestoShowListBoxMargFSPvT.Items{iSpecies})];
    end
end

speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);

% Compute the marginal distributions at all times and for all species.
Nd = app.SSITModel.Solutions.fsp{1}.p.dim;
mdist = cell(length(T_array2),nSpecies);
for it = 1:length(T_array2)
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


if ~app.FspMarginalTimeCreateMovieCheckBox.Value
    for iplt = species2Plot
        figure()
        for it = 1:length(T_array2)
            stairs([0:length(mdist{it,iplt})], [mdist{it,iplt};0],'linewidth',2);
            hold('on');
        end
        title(['Marginals of ',char(app.SpeciestoShowListBoxMargFSPvT.Items(iplt))])
        xlabel('Species Count')
        ylabel('Probability')
        legendCell = cellstr(num2str(T_array2', 'Time=%1.2f'));
        legend(legendCell)
        set(gca,'fontsize',20)
    end
else
    fig_save = figure();
    str = ['TimeMarginal'];
    v = VideoWriter(str,'MPEG-4');
    open(v)
    Nd = app.SSITModel.Solutions.fsp{1}.p.dim;
    for it = 1:length(T_array2)
        hold off
        for iplt = species2Plot
            stairs([0:length(mdist{it,iplt})], [mdist{it,iplt};0],'linewidth',2);
            hold on

        end
        title(['Marginals Distributions'])
        xlabel('Species Count')
        ylabel('Probability')
        set(gca,'fontsize',20)
        YLIM = get(gca,'ylim');YLIM(1)=0;
        set(gca,'ylim',YLIM);

        F1 = getframe(gcf);
        writeVideo(v,F1)
    end
    close(v)
    save_button = uicontrol;
    save_button.String = 'Save Movie';
    save_button.Callback = eval(['@SaveMovie']);
end

    function SaveMovie(~,~)
        Savedir = fullfile(pwd,['TimeMarginal.mp4']);
        [name,pathname] = uiputfile('*.mp4', 'Save movie as');
        MovieName = [pathname,name];
        movefile(Savedir,MovieName)
    end
end
