function loadDataFromFile(app,fname,pathname)
%% This code will prompt the user to choose a file to load data from disk. 
% Then, it will load the data to memory.
% The code will find the appropriate column headers and use these to
%      populate the dropdown menus to connect data to species in the GUI.

if nargin<3
    [fname,pathname] = uigetfile('*.xlsx;*.xls;*.csv','Pick the spreadsheet with your data');
end

app.DataLoadingAndFittingTabOutputs.dataFileName = [pathname,fname];
app.DataFileNoneSelectedLabel.Text = ['Data File: ',pathname,fname];

Tab = readtable([pathname,fname]);

app.TotalCellsInDataLabel.Text = ['Total Cells in Data: ',num2str(size(Tab,1))];

columns = Tab.Properties.VariableNames;
app.DataLoadingAndFittingTabOutputs.columns = columns;
txt = {};
app.DataLoadingAndFittingTabOutputs.columnTypes = {};
for iField = 1:length(columns)
    if isnumeric(Tab{1,iField})
        txt{end+1} = [columns{iField},': ',num2str(Tab{1,iField}),', numeric'];
        app.DataLoadingAndFittingTabOutputs.columnTypes{iField} = 'numeric';
    elseif iscell(Tab{1,iField})
        txt{end+1} = [columns{iField},': ',Tab{1,iField}{1},', string'];
        app.DataLoadingAndFittingTabOutputs.columnTypes{iField} = 'string';
    elseif ischar(Tab{1,iField})
        txt{end+1} = [columns{iField},': ',Tab{1,iField},', string'];
        app.DataLoadingAndFittingTabOutputs.columnTypes{iField} = 'string';
    end        
end
app.FieldsinDataTextArea.Value = txt;

J = strcmp(columns,'time');
if sum(J) == 0
    options = columns;
    [idx, tf] = listdlg('PromptString', 'Select the column for time:', ...
        'SelectionMode', 'single', ...
        'ListString', options);
    Tab.Properties.VariableNames{idx} = 'time';
end

app.ParEstFitTimesList.Items = cellstr(num2str(unique(Tab.time)));

nSp = length(app.SSITModel.species);
for iSp = 1:nSp
    app.(['DataSpecies',num2str(iSp)]).Items = ['none';columns'];
end

for iConstr = 1:6
    if iConstr>1
        app.(['AndOr',num2str(iConstr)]).Value = '-';
    end
    app.(['DataConstrText',num2str(iConstr)]).Value = '-';

    app.(['DataConstrChoice',num2str(iConstr)]).Items = ['-';columns'];
    app.(['DataConstrChoice',num2str(iConstr)]).Value = '-';
end

% for 
% DATA = table2cell(Tab);
% app.DataLoadingAndFittingTabOutputs.columns = TXT;
% app.DataLoadingAndFittingTabOutputs.dataTable = DATA;
% 
% Q = contains(TXT,'time')|contains(TXT,'Time')|contains(TXT,'TIME');
% 
% if sum(Q)==1
%     app.ParEstFitTimesList.Items = {};
%     app.ParEstFitTimesList.Value = {};
%     col_time = find(Q);
%     app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_time_index = col_time;
%     app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times = sort(unique(cell2mat(DATA(:,col_time))));
%     for i=1:length(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)
%         app.ParEstFitTimesList.Items{i} = num2str(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times(i));
%         app.ParEstFitTimesList.Value{i} = num2str(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times(i));
%     end
%     maxFitTime = max(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times);
%     fitTimes = sort(unique([(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)',(str2num(app.FspPrintTimesField.Value))]));
%     fitTimes = fitTimes(fitTimes<=maxFitTime);
%     app.FspPrintTimesField.Value = ['[',num2str(fitTimes),']'];
%     % We need to make sure that the fitting times are included in the solution times.
% 
%     app.ParEstFitTimesList.Value = app.ParEstFitTimesList.Items;  % default is to select all times.
% 
% end
% 
% app.DataLoadingAndFittingTabOutputs.dataTable = DATA;
end

