function makeTimeMarginalsSensTab(app)
% This function creates plots comparing the different marginals throughout
% time of a single species. ReactionsTabOutputs.inputs require vector format for the time
% points.
T_array = eval(app.SensPrintTimesEditField.Value);
T_array2 = eval(app.SensMarginalTimeVectorEditField_2.Value);
Plts_to_make = [app.SensMarginalTimeX1CheckBox.Value,app.SensMarginalTimeX2CheckBox.Value,app.SensMarginalTimeX3CheckBox.Value];
solutionFormat = app.SensFspTabOutputs.solutions.format;
if app.SensMarginalTimeCreateMovieCheckBox.Value == 0
    for iplt = 1:3
        if Plts_to_make(iplt)
            figure()
            for i = 1:length(T_array2)
                [~,j] =  min(abs(T_array-T_array2(i)));
                % Compute the marginal distributions
                if (strcmp(solutionFormat, 'forward'))
                    mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states, ...
                        app.SensFspTabOutputs.solutions.data{j}.p,false);
                else
                    mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.ps{j}.states, ...
                        app.SensFspTabOutputs.solutions.data.ps{j}.p,false);
                end
                stairs([0:length(mdist{iplt})], [mdist{iplt};0],'linewidth',2);
                hold('on');
            end
            title(['Marginals of ',char(app.NameTable.Data(iplt,2))])
            xlabel('Species Count')
            ylabel('Probability')
            legendCell = cellstr(num2str(T_array2', 'Time=%d'));
            legend(legendCell)
            set(gca,'fontsize',20)
        end
    end
else
    for iplt = 1:3
        if Plts_to_make(iplt)
            str = ['DefaultTimeMarginalX',num2str(iplt)];
            v{iplt} = VideoWriter(str);
            open(v{iplt})
            fig_save = figure();
            for i = 1:length(T_array2)
                [~,j] =  min(abs(T_array-T_array2(i)));
                if (strcmp(solutionFormat, 'forward'))
                    mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states, ...
                        app.SensFspTabOutputs.solutions.data{j}.p,false);
                else
                    mdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.ps{j}.states, ...
                        app.SensFspTabOutputs.solutions.data.ps{j}.p,false);
                end
                stairs([0:length(mdist{iplt})], [mdist{iplt};0],'linewidth',2);
                title(['Marginals of x',num2str(iplt)])
                xlabel('Species Count')
                ylabel('Probability')
                set(gca,'fontsize',20)

                F1 = getframe(gcf);
                writeVideo(v{iplt},F1)
            end
            close(v{iplt})
            save_button = uicontrol;
            save_button.String = 'Save Movie';
            save_button.Callback = eval(['@SaveMovieX',num2str(iplt)]);
        end
    end
end
    function SaveMovieX1(~,~)
        Savedir = fullfile(pwd,'DefaultTimeMarginalX1.avi');
        prompt = 'File Name? (End in .avi)';
        opts.Resize = 'on';
        definput = {'DefaultTimeMarginalX1.avi'};
        name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
        MovieName = string(name);
        movefile(Savedir,MovieName)
    end
    function SaveMovieX2(~,~)
        Savedir = fullfile(pwd,'DefaultTimeMarginalX2.avi');
        prompt = 'File Name? (End in .avi)';
        opts.Resize = 'on';
        definput = {'DefaultTimeMarginalX2.avi'};
        name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
        MovieName = string(name);
        movefile(Savedir,MovieName)
    end
    function SaveMovieX3(~,~)
        Savedir = fullfile(pwd,'DefaultTimeMarginalX3.avi');
        prompt = 'File Name? (End in .avi)';
        opts.Resize = 'on';
        definput = {'DefaultTimeMarginalX3.avi'};
        name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
        MovieName = string(name);
        movefile(Savedir,MovieName)
    end
end
