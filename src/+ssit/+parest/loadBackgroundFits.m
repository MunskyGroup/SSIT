function loadBackgroundFits(app,choice)
%% Put results of the fits back into the app.
load([app.DataLoadingAndFittingTabOutputs.backgroundFitDir,'/',choice],'Loc_app_str','fields');

for i=1:length(fields)
    try
        eval(['app.',fields{i},' = Loc_app_str.',fields{i},';'])
    catch
    end
end

ssit.parest.updateModelSolveAndCompareToData(app);
makeSeparatePlotOfData(app);
figure(app.UIFigure);
