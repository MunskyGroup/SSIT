function loadDataFromFile(app,fname,pathname)
%% This code will prompt the user to choose a file to load data from disk. 
% Then, it will load the data to memory.
% The code will find the appropriate column headers and use these to
%      populate the dropdown menus to connect data to species in the GUI.

if nargin<3
    [fname,pathname] = uigetfile('*.xlsx;*.xls;*.csv','Pick the spreadsheet with your data');
end

Tab = readtable([pathname,fname]);
TXT = Tab.Properties.VariableNames;
DATA = table2cell(Tab);

app.DataLoadingAndFittingTabOutputs.columns = TXT;
app.DataLoadingAndFittingTabOutputs.dataTable = DATA;

Q = contains(TXT,'time')|contains(TXT,'Time')|contains(TXT,'TIME');

if sum(Q)==1
    app.ParEstFitTimesList.Items = {};
    app.ParEstFitTimesList.Value = {};
    col_time = find(Q);
    app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_time_index = col_time;
    app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times = sort(unique(cell2mat(DATA(:,col_time))));
    for i=1:length(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)
        app.ParEstFitTimesList.Items{i} = num2str(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times(i));
        app.ParEstFitTimesList.Value{i} = num2str(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times(i));
    end
    maxFitTime = max(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times);
    fitTimes = sort(unique([(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)',(str2num(app.FspPrintTimesField.Value))]));
    fitTimes = fitTimes(fitTimes<=maxFitTime);
    app.FspPrintTimesField.Value = ['[',num2str(fitTimes),']'];
    % We need to make sure that the fitting times are included in the solution times.
   
    app.ParEstFitTimesList.Value = app.ParEstFitTimesList.Items;  % default is to select all times.

end

app.DataLoadingAndFittingTabOutputs.dataTable = DATA;
end

