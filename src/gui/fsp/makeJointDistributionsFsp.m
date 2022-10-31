function makeJointDistributionsFsp(app)
% This function creates the joint distributions of the FSP results as
% either a movie or not. Also either provides a contour plot or a mesh plot
% depending on the selection of the user.
nD = app.FspTabOutputs.solutions{1}.p.dim;
if nD<2
    error('Only one species avaialble -- cannot make joint pdf')
end
T_array = eval(app.FspPrintTimesField.Value);

allInds = [1:nD];
sp1 = find(strcmp(app.ReactionsTabOutputs.varNames,app.JointSp1.Value));
sp2 = find(strcmp(app.ReactionsTabOutputs.varNames,app.JointSp2.Value));

if nD==2
    for it = length(T_array):-1:1
        Joints{it} = double(app.FspTabOutputs.solutions{it}.p.data);
    end
else
    for it = length(T_array):-1:1
        k = setdiff(allInds,[sp1,sp2]);
        Joints{it} = double(app.FspTabOutputs.solutions{it}.p.sumOver(k).data);
    end
end

titleTxt = ['Joint distribution of ',app.ReactionsTabOutputs.varNames{sp1},' and ',app.ReactionsTabOutputs.varNames{sp2}];
yLabTxt = ['Species ',app.ReactionsTabOutputs.varNames{sp1}];
xLabTxt = ['Species ',app.ReactionsTabOutputs.varNames{sp2}];
if app.FspJointCreateMovCheckBox.Value == 0    % Creates non-movie plot at the selected timepoint
    [~,j] = min(abs(T_array-app.FspTimeSlider.Value));
    figure()
    Z = log10((Joints{j}));
    switch app.FspContourCheckBox.Value
        case 0
            mesh(0:size(Z,2)-1,0:size(Z,1)-1,Z)
        case 1
            contourf(0:size(Z,2)-1,0:size(Z,1)-1,Z)
            colorbar
    end
    title(titleTxt);ylabel(yLabTxt);xlabel(xLabTxt);
    set(gca,'fontsize',20);
elseif app.FspJointCreateMovCheckBox.Value == 1    % Creates non-movie plot at the selected timepoint
    figure()
    str = 'DefaultJoint';
    v = VideoWriter(str,'MPEG-4');
    open(v)
    for j = 1:length(T_array)
        Z = log10((Joints{j}));
        switch app.FspContourCheckBox.Value
            case 0
                mesh(0:size(Z,2)-1,0:size(Z,1)-1,Z)
            case 1
                contourf(0:size(Z,2)-1,0:size(Z,1)-1,Z)
                colorbar
        end
        title(titleTxt);ylabel(yLabTxt);xlabel(xLabTxt)
        set(gca,'fontsize',20);
        F1 = getframe(gcf);
        writeVideo(v,F1)
    end
    close(v)
    save_button = uicontrol;
    save_button.String = 'Save Movie';
    save_button.Callback = @(a,b)SaveMovie(a,b);

end
end

function SaveMovie(~,~)
Savedir = fullfile(pwd,'DefaultJoint.mp4');
[name,pathname] = uiputfile('*.mp4', 'Save movie as');
MovieName = [pathname,name];
movefile(Savedir,MovieName)
end