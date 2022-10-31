function makePlotOfData(app)
%% This function creates a histogram from the loaded data.
error('OLD CODE - PLEASE REMOVE CALL')
return
it = find(strcmp(app.ParEstPlotTimesDropDown.Items,app.ParEstPlotTimesDropDown.Value));
hold(app.data_histogram_plot,'off')
L = {};
xm = 0; ym = 0;

% Switches plotted data, (x1, x2, x3), and legend depending on if the
% check box was selected
for icb=1:3
    switch icb
        case 1
            cb = app.ParEstX1CheckBox.Value;
            so = [app.DataLoadingAndFittingTabOutputs.xInd(2) app.DataLoadingAndFittingTabOutputs.xInd(3)];
            le = app.NameTable.Data{1,2};
        case 2
            cb = app.ParEstX2CheckBox.Value;
            so = [app.DataLoadingAndFittingTabOutputs.xInd(1) app.DataLoadingAndFittingTabOutputs.xInd(3)];
            le = app.NameTable.Data{1,2};
        case 3
            cb = app.ParEstX3CheckBox.Value;
            so = [app.DataLoadingAndFittingTabOutputs.xInd(1) app.DataLoadingAndFittingTabOutputs.xInd(2)];
            le = app.NameTable.Data{1,2};
    end
    
    if cb
        matTensor = double(app.DataLoadingAndFittingTabOutputs.dataTensor);
        H1 = matTensor(it,:);
        H1(end+1) = 0;
        H1 = H1/sum(H1);  % Normalize data before plotting.
        
        if app.DensityButton_2.Value
            stairs(app.data_histogram_plot,[0:length(H1)-1],H1);
            ym = max(ym,max(H1));
        else
            stairs(app.data_histogram_plot,[0:length(H1)-1],cumsum(H1));
            ym=1;
        end
        hold(app.data_histogram_plot,'on')
        L{end+1} = [le,'-data'];
        xm = max(xm,length(H1));
        
        try
            so = so(so~=0);
            if ~isempty(so)
                H1 = squeeze(sum(app.DataLoadingAndFittingTabOutputs.fitResults.current(it,:,:,:),so));
            else
                H1 = squeeze(app.DataLoadingAndFittingTabOutputs.fitResults.current(it,:,:,:));
            end
            H1(end+1)=0;
            if app.DensityButton_2.Value
                stairs(app.data_histogram_plot,[0:length(H1)-1],H1,'linewidth',3);
                ym = max(ym,max(H1));
            else
                stairs(app.data_histogram_plot,[0:length(H1)-1],cumsum(H1),'linewidth',3);
                ym=1;
            end
            hold(app.data_histogram_plot,'on')
            L{end+1} = [le,'-mod'];
        catch
        end
        
    end
end


%% Make trajectory plot for model.
T_array = eval(app.FspPrintTimesField.Value);
for it = 1:length(T_array)
    if ~isempty(app.FspTabOutputs.solutions{it})
        mdist = ssit.fsp.marginals(app.FspTabOutputs.solutions{it}.states, app.FspTabOutputs.solutions{it}.p);
        for j=1:3
            mns(it,j) = [0:length(mdist{j})-1]*mdist{j};
            mns2(it,j) = [0:length(mdist{j})-1].^2*mdist{j};
        end
    end
end
vars = mns2-mns.^2;
Plts_to_make = [app.ParEstX1CheckBox.Value,app.ParEstX2CheckBox.Value,app.ParEstX3CheckBox.Value];
cols = ['b','r','g'];
cols2 = [.90 .90  1.00; 1.00 .90 .90; .90 1.00 .90];
LG = {};
hold(app.data_traj_plot,'off');
for iplt=1:3
    if Plts_to_make(iplt)
        BD = [mns(:,iplt)'+sqrt(vars(:,iplt)'),mns(end:-1:1,iplt)'-sqrt(vars(end:-1:1,iplt)')];
        TT = [T_array(1:end),T_array(end:-1:1)];
        fill(app.data_traj_plot,TT,BD,cols2(iplt,:));
        hold(app.data_traj_plot,'on');
        plot(app.data_traj_plot,T_array,mns(:,iplt),cols(iplt),'linewidth',2);
        LG{end+1} = [char(app.NameTable.Data(iplt,2)),' Model Mean \pm std'];
        LG{end+1} = [char(app.NameTable.Data(iplt,2)),' Model Mean'];
    end
end
app.data_traj_plot.XLim = [min(T_array),max(T_array)];

%% Add data to trajectory plot.
for iplt = 1:3
    switch iplt
        case 1
            so = [app.DataLoadingAndFittingTabOutputs.xInd(2) app.DataLoadingAndFittingTabOutputs.xInd(3)];
        case 2
            so = [app.DataLoadingAndFittingTabOutputs.xInd(1) app.DataLoadingAndFittingTabOutputs.xInd(3)];
        case 3
            so = [app.DataLoadingAndFittingTabOutputs.xInd(1) app.DataLoadingAndFittingTabOutputs.xInd(2)];
    end
    if Plts_to_make(iplt)
        matTensor = double(app.DataLoadingAndFittingTabOutputs.dataTensor);
        P = matTensor./sum(matTensor,2);
        mn = P*[0:size(matTensor,2)-1]';
        mn2 = P*([0:size(matTensor,2)-1]').^2;
        var = mn2-mn.^2;
        plot(app.data_traj_plot,app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,mn,[cols(iplt),...
            'o'],'MarkerSize',12,'MarkerFaceColor',cols(iplt))
        errorbar(app.data_traj_plot,app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,...
            mn,sqrt(var),[cols(iplt),...
            'o'],'MarkerSize',12,'MarkerFaceColor',cols(iplt),'linewidth',3)
        LG{end+1} = [char(app.NameTable.Data(iplt,2)),' Data'];
    end
end
legend(app.data_traj_plot,LG)
        %%

try
    app.data_histogram_plot.XTick = round(linspace(0,xm,5),2,'significant');
    app.data_histogram_plot.XLim = [0,xm];
    app.data_histogram_plot.YLim = [0,ym];
    legend(app.data_histogram_plot,L);
catch
end
end