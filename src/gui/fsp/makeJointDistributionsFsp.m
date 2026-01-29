function makeJointDistributionsFsp(app)
% This function creates the joint distributions of the FSP results as
% either a movie or not. Also either provides a contour plot or a mesh plot
% depending on the selection of the user.
nD = app.SSITModel.Solutions.fsp{1}.p.dim;
if nD<2
    error('Only one species avaialble -- cannot make joint pdf')
end
T_array = eval(app.FspPrintTimesField.Value);

sp1 = find(strcmp(app.JointSp1.Items,app.JointSp1.Value));
sp2 = find(strcmp(app.JointSp1.Items,app.JointSp2.Value));

% Check if plot selection is distorted.
if sp1>nD
    sp1p = find(strcmp(app.JointSp1.Items,app.JointSp1.Value(1:end-12)));
    sp1pp = sp1 - nD;
    distort1 = true;
else
    sp1p = sp1;
    sp1pp = sp1;
    distort1 = false;
end
if sp2>nD
    sp2p = find(strcmp(app.JointSp2.Items,app.JointSp2.Value(1:end-12)));
    sp2pp = sp2-nD;
    distort2 = true;
else
    sp2p = sp2;
    sp2pp = sp2;
    distort2 = false;
end

% if sp1p==sp2p
%     % error('Only one species chosen -- cannot make joint pdf')
% end

% if nD==2
%     for it = length(T_array):-1:1
%         Joints{it} = double(app.SSITModel.Solutions.fsp{it}.p.data);
%     end
% else
allInds = [1:nD];
for it = length(T_array):-1:1
    k = setdiff(allInds,[sp1p,sp2p]);
    Joints{it} = double(app.SSITModel.Solutions.fsp{it}.p.sumOver(k).data);
    if sp1p==sp2p
        Joints{it} = diag(Joints{it});
    end
    if distort2
        Joints{it} = app.SSITModel.pdoOptions.PDO.conditionalPmfs{sp2pp}(:,1:length(Joints{it}))*Joints{it};
    end
    if distort1
        Joints{it} = Joints{it}*(app.SSITModel.pdoOptions.PDO.conditionalPmfs{sp1pp}(:,1:length(Joints{it})))';
    end
end
% end

titleTxt = ['Joint distribution of ',app.JointSp1.Items{sp1},' and ',app.JointSp1.Items{sp2}];
yLabTxt = [app.JointSp1.Items{sp1}];
xLabTxt = [app.JointSp1.Items{sp2}];
if app.FspJointCreateMovCheckBox.Value == 0    % Creates non-movie plot at the selected timepoint
    [~,j] = min(abs(T_array-app.FspTimeSlider.Value));
    figure()
    Z = log10(max(1e-8,Joints{j}));
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
        Z = log10(max(1e-8,Joints{j}));
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