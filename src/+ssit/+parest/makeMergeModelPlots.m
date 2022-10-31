function [FitResults] = makeMergeModelPlots(app,SSITapp,callFromModelMerge)
arguments
    app
    SSITapp
    callFromModelMerge = false
end

vars_to_fit = strcmp(app.fit_parameters_table.Data(:,3),'y');
x0 = log10(cell2mat(app.fit_parameters_table.Data(vars_to_fit,2)));
FitResults = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FitResults;

k1 = find(strcmp(app.MergedModelDropDown.Items, app.MergedModelDropDown.Value));
k2 = find(strcmp(app.MergedModelDropDown2.Items, app.MergedModelDropDown2.Value));

if isempty(app.DataLoadingAndFittingTabOutputs.mergeFitResults)
    for i=1:length(FitResults)
        app.DataLoadingAndFittingTabOutputs.mergeFitResults{i} = FitResults{i}(x0);
        fspBounds{i} = app.DataLoadingAndFittingTabOutputs.mergeFitResults{i}{4};
    end
end

% update choices for time based on user selection.
for j = 1:length(app.DataLoadingAndFittingTabOutputs.mergeFitResults{k1}{3})
    app.ParEstPlotTimesDropDown.Items{j} =num2str(app.DataLoadingAndFittingTabOutputs.mergeFitResults{k1}{3}(j));
end
app.ParEstPlotTimesDropDown.Items=app.ParEstPlotTimesDropDown.Items(1:j);

for j = 1:length(app.DataLoadingAndFittingTabOutputs.mergeFitResults{k2}{3})
    app.ParEstPlotTimesDropDown_2.Items{j} =num2str(app.DataLoadingAndFittingTabOutputs.mergeFitResults{k2}{3}(j));
end
app.ParEstPlotTimesDropDown_2.Items=app.ParEstPlotTimesDropDown_2.Items(1:j);


% Make plot for chosen model and time point.
it1 = find(strcmp(app.ParEstPlotTimesDropDown.Items,app.ParEstPlotTimesDropDown.Value));
Nmax = size(app.DataLoadingAndFittingTabOutputs.mergeFitResults{k1}{1},2)-1;
hold(app.data_histogram_plot,'off')
stairs(app.data_histogram_plot,[0:Nmax],app.DataLoadingAndFittingTabOutputs.mergeFitResults{k1}{1}(it1,:),'r','linewidth',2)
hold(app.data_histogram_plot,'on')
stairs(app.data_histogram_plot,[0:Nmax],app.DataLoadingAndFittingTabOutputs.mergeFitResults{k1}{2}(it1,:),'b','linewidth',2)

% Make plot for chosen model and time point.
it2 = find(strcmp(app.ParEstPlotTimesDropDown_2.Items,app.ParEstPlotTimesDropDown_2.Value));
Nmax = size(app.DataLoadingAndFittingTabOutputs.mergeFitResults{k2}{1},2)-1;
hold(app.data_histogram_plot2,'off')
stairs(app.data_histogram_plot2,[0:Nmax],app.DataLoadingAndFittingTabOutputs.mergeFitResults{k2}{1}(it2,:),'r','linewidth',2)
hold(app.data_histogram_plot2,'on')
stairs(app.data_histogram_plot2,[0:Nmax],app.DataLoadingAndFittingTabOutputs.mergeFitResults{k2}{2}(it2,:),'b','linewidth',2)

% Reset Merege model with new FSP constraints.
for i=1:length(FitResults)
    fspBounds{i} = app.DataLoadingAndFittingTabOutputs.mergeFitResults{i}{4};
end

if ~callFromModelMerge
    ssit.parest.ModelMergeMerge(app,SSITapp,false,fspBounds,false);
end

end