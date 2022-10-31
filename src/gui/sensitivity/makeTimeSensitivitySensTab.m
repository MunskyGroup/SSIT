function makeTimeSensitivitySensTab(app)
% This function creates plots comparing the different marginals throughout
% time of a single species. ReactionsTabOutputs.inputs require vector format for the time
% points.
T_array = eval(app.SensPrintTimesEditField.Value);
T_array2 = eval(app.SensMarginalTimeVectorEditField_2.Value);
Plts_to_make = [app.SensMarginalTimeX1CheckBox.Value,app.SensMarginalTimeX2CheckBox.Value,app.SensMarginalTimeX3CheckBox.Value];
% Find the parameter index to plot
ipar = 1;
num_pars = length(app.ReactionsTabOutputs.parameters);
while (~strcmp(app.ReactionsTabOutputs.parameters{ipar}, app.SensParDropDown.Value) && ipar < num_pars)
    ipar = ipar + 1;
end
solutionFormat = app.SensFspTabOutputs.solutions.format;
if app.SensMarginalTimeCreateMovieCheckBox.Value == 0
    for iplt = 1:3
        if Plts_to_make(iplt)
            figure()
            for i = 1:length(T_array2)
                [~,j] =  min(abs(T_array-T_array2(i)));
                % Compute the Sensitivity Distributions
                if (strcmp(solutionFormat, 'forward'))
                    sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states,...
                        app.SensFspTabOutputs.solutions.data{j}.dp(:, ipar),false);
                else
                    sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.dps{j, ipar}.states,...
                        app.SensFspTabOutputs.solutions.data.dps{j, ipar}.p,false);
                end
                stairs([0:length(sensmdist{iplt})], [sensmdist{iplt};0],'linewidth',2);
                hold('on');
            end
            paramTitle = horzcat(' w.r.t. ',app.SensParDropDown.Value);
            strTitle = horzcat('Sensitivity of x',num2str(iplt));
            titleTitle = append(strTitle,paramTitle);
            title(titleTitle)
            xlabel('Species Count')
            ylabel('Sensitivity')
            legendCell = cellstr(num2str(T_array2', 'Time=%d'));
            legend(legendCell)
            set(gca,'fontsize',18)
        end
    end
else
    for iplt = 1:3
        if Plts_to_make(iplt)
            str = ['DefaultTimeSensitivityX',num2str(iplt)];
            v{iplt} = VideoWriter(str);
            open(v{iplt})
            fig_save = figure();
            for i = 1:length(T_array2)
                [~,j] =  min(abs(T_array-T_array2(i)));
                %Compute the Sensitivity Distributions
                if (strcmp(solutionFormat, 'forward'))
                    sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{j}.states,...
                        app.SensFspTabOutputs.solutions.data{j}.dp(:, ipar),false);
                else
                    sensmdist = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.dps{j, ipar}.states,...
                        app.SensFspTabOutputs.solutions.data.dps{j, ipar}.p,false);
                end
                stairs([0:length(sensmdist{iplt})], [sensmdist{iplt};0],'linewidth',2);
                paramTitle = horzcat(' w.r.t. ',app.SensParDropDown.Value);
                strTitle = horzcat('Sensitivity of x',num2str(iplt));
                titleTitle = append(strTitle,paramTitle);
                title(titleTitle)
                xlabel('Species Count')
                ylabel('Sensitivity')
                set(gca,'fontsize',18)
                
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
        Savedir = fullfile(pwd,'DefaultTimeSensitivityX1.avi');
        prompt = 'File Name? (End in .avi)';
        opts.Resize = 'on';
        definput = {'DefaultTimeSensitivityX1.avi'};
        name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
        MovieName = string(name);
        movefile(Savedir,MovieName)
    end
    function SaveMovieX2(~,~)
        Savedir = fullfile(pwd,'DefaultTimeSensitivityX2.avi');
        prompt = 'File Name? (End in .avi)';
        opts.Resize = 'on';
        definput = {'DefaultTimeSensitivityX2.avi'};
        name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
        MovieName = string(name);
        movefile(Savedir,MovieName)
    end
    function SaveMovieX3(~,~)
        Savedir = fullfile(pwd,'DefaultTimeSensitivityX3.avi');
        prompt = 'File Name? (End in .avi)';
        opts.Resize = 'on';
        definput = {'DefaultTimeSensitivityX3.avi'};
        name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
        MovieName = string(name);
        movefile(Savedir,MovieName)
    end
end