function makeJointDistributionsSens(app)
% This function creates the joint distributions of the FSP results as
% either a movie or not. Also either provides a contour plot or a mesh plot
% depending on the selection of the user.

T_array = eval(app.SensPrintTimesEditField.Value);
solutionFormat = app.SensFspTabOutputs.solutions.format;

for it = length(T_array):-1:1
    % Compute the marginal distributions
    if (strcmp(solutionFormat, 'forward'))
        [~,joints] = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data{it}.states, ...
            app.SensFspTabOutputs.solutions.data{it}.p,false);
    else
        [~,joints] = ssit.fsp.marginals(app.SensFspTabOutputs.solutions.data.ps{it}.states, ...
            app.SensFspTabOutputs.solutions.data.ps{it}.p,false);
    end
    Joints{it} = joints;
end

plts = [app.SensX1X2CheckBox.Value,app.SensX1X3CheckBox.Value,app.SensX2X3CheckBox.Value];
pltsinds = [1 2;1 3;2 3];

if app.SensCreateMovieCheckBox.Value == 0    % Creates non-movie plot at the selected timepoint
    [~,j] = min(abs(T_array-app.SensPlotTimeSlider.Value));
    for iplt=1:length(plts)
        if plts(iplt)
            figure()
            Z = log10((Joints{j}{pltsinds(iplt,1),pltsinds(iplt,2)}));
            switch app.SensContourCheckBox.Value
                case 0
                    mesh(Z)
                case 1
                    contourf(Z)
                    colorbar
            end
            title(['Joint distribution of x',num2str(pltsinds(iplt,1)),' and x',num2str(pltsinds(iplt,2))]);
            ylabel(['Species x',num2str(pltsinds(iplt,1))])
            xlabel(['Species x',num2str(pltsinds(iplt,2))])
            set(gca,'fontsize',20);
        end
    end
elseif app.SensCreateMovieCheckBox.Value == 1    % Creates non-movie plot at the selected timepoint
    for iplt=1:length(plts)
        if plts(iplt)
            figure()
            str = ['DefaultJointX',num2str(pltsinds(iplt,1)),'X',num2str(pltsinds(iplt,2))];
            v{iplt} = VideoWriter(str);
            open(v{iplt})
            for j = 1:length(T_array)
                Z = log10((Joints{j}{pltsinds(iplt,1),pltsinds(iplt,2)}));
                switch app.SensContourCheckBox.Value
                    case 0
                        mesh(Z)
                    case 1
                        contourf(Z)
                        colorbar
                end
                title(['Joint distribution of x',num2str(pltsinds(iplt,1)),' and x',num2str(pltsinds(iplt,2))]);
                ylabel(['Species x',num2str(pltsinds(iplt,1))])
                xlabel(['Species x',num2str(pltsinds(iplt,2))])
                set(gca,'fontsize',20);
                F1 = getframe(gcf);
                writeVideo(v{iplt},F1)
            end
            close(v{iplt})
            save_button = uicontrol;
            save_button.String = 'Save Movie';
            save_button.Callback = eval(['@SaveMovie',num2str(iplt)]);
        end
    end
end

end

function SaveMovie1(~,~)
Savedir = fullfile(pwd,'DefaultJointX1X2.avi');
prompt = 'File Name? (End in .avi)';
opts.Resize = 'on';
definput = {'DefaultJointX1X2.avi'};
name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
MovieName = string(name);
movefile(Savedir,MovieName)
end
function SaveMovie2(~,~)
Savedir = fullfile(pwd,'DefaultJointX1X3.avi');
prompt = 'File Name? (End in .avi)';
opts.Resize = 'on';
definput = {'DefaultJointX1X3.avi'};
name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
MovieName = string(name);
movefile(Savedir,MovieName)
end
function SaveMovie3(~,~)
Savedir = fullfile(pwd,'DefaultJointX2X3.avi');
prompt = 'File Name? (End in .avi)';
opts.Resize = 'on';
definput = {'DefaultJointX2X3.avi'};
name = inputdlg(prompt,'Save Movie',[1 40],definput,opts);
MovieName = string(name);
movefile(Savedir,MovieName)
end
