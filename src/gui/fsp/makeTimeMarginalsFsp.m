function makeTimeMarginalsFsp(app)
% This function creates plots comparing the different marginals throughout
% time of a single species. ReactionsTabOutputs.inputs require vector format for the time
% points.

T_array = eval(app.FspPrintTimesField.Value);
T_array2 = eval(app.FspMarginalVecField.Value);
Nd = app.SSITModel.Solutions.fsp{1}.p.dim;
for iSpecies = 1:Nd
    Plts_to_make(iSpecies) = max(contains(app.SpeciestoShowListBoxMargFSPvT.Value,app.SSITModel.species{iSpecies}));
end
if ~app.FspMarginalTimeCreateMovieCheckBox.Value
    for iplt = 1:Nd
        if Plts_to_make(iplt)
            INDS = setdiff([1:Nd],iplt);
            figure()
            for i = 1:length(T_array2)
                [~,j] =  min(abs(T_array-T_array2(i)));
                % Compute the marginal distributions
                if Nd==1
                    mdist = double(app.SSITModel.Solutions.fsp{j}.p.data);
                else
                    mdist = double(app.SSITModel.Solutions.fsp{j}.p.sumOver(INDS).data);
                end
                stairs([0:length(mdist)], [mdist;0],'linewidth',2);
                hold('on');
            end
            title(['Marginals of ',char(app.SSITModel.species(iplt))])
            xlabel('Species Count')
            ylabel('Probability')
            legendCell = cellstr(num2str(T_array2', 'Time=%1.2f'));
            legend(legendCell)
            set(gca,'fontsize',20)
        end
    end
else
    fig_save = figure();
    str = ['TimeMarginal'];
    v = VideoWriter(str,'MPEG-4');
    open(v)
    Nd = app.SSITModel.Solutions.fsp{1}.p.dim;
    for i = 1:length(T_array2)
        hold off
        [~,j] =  min(abs(T_array-T_array2(i)));
        for iplt = 1:Nd
            if Plts_to_make(iplt)
                if Nd==1
                    mdist{1} = double(app.SSITModel.Solutions.fsp{j}.p.data);
                else
                    for ii=1:Nd
                        INDS = setdiff([1:Nd],ii);
                        mdist{ii} = double(app.SSITModel.Solutions.fsp{j}.p.sumOver(INDS).data);
                    end
                end
                stairs([0:length(mdist{iplt})], [mdist{iplt};0],'linewidth',2);
                hold on
            end
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
