function marginalizeDataGui(app)
% This functions opens an interface so the user can match columns in their
% data file to species in their models.

% Store data from previous app
colHeaders = app.DataLoadingAndFittingTabOutputs.columns;
xcelData = app.DataLoadingAndFittingTabOutputs.dataTable;
xcelData = xcelData([1:1:10],:);
Nd = size(app.NameTable.Data,1);
rowHeaders = {app.NameTable.Data{:,2},'time','Counts/Probability','Marginalize over','Condition on ='};

app.DataLoadingAndFittingTabOutputs.marginalMatrix = zeros(length(rowHeaders)-1,length(colHeaders));
app.DataLoadingAndFittingTabOutputs.conditionOnArray = zeros(1,length(colHeaders));
% Create the figure
hFig = figure('Position',[0 0 1300 800],'Visible','off');
movegui(hFig,'center')
set(hFig,'Tag','marginalizeDataGui','HandleVisibility','on','Visible','on');
% Create Table Snapshot
ssTable = uitable(hFig);
ssTable.Data = xcelData;
ssTable.ColumnName = colHeaders;
ssTable.Position = [40 650 700 112];
ssTable.FontSize = 16;
uicontrol(hFig,'Style','text','String','Snapshot of Loaded Data',...
        'position',[60 770 400 22],'fontsize',16)

% Create labels
for iH = 1:length(colHeaders)
    uicontrol(hFig,'Style','text','String',colHeaders{iH},...
        'position',[50+140*(iH) 610 100 22],'fontsize',16)
end
for iR = 1:length(rowHeaders)
    uicontrol(hFig,'Style','text','String',rowHeaders{iR},...
        'position',[20 570-(40*(iR-1)) 170 22],'fontsize',16)
end

% Create Checkboxes
for kR = length(rowHeaders) - 1 : -1 : 1
    for kH = length(colHeaders) : -1 :1
        cbh(kR,kH) = uicontrol(hFig,'Style','checkbox', ...
            'Value',0,'Position',[50+140*kH+40 570-(40*(kR-1)) 100 22]);
    end
end

% auto-detect and click 'time'
Itime = Nd+1; %find(contains(rowHeaders,'time')|contains(rowHeaders,'Time'));
Jtime = find(contains(colHeaders,'time')|contains(colHeaders,'Time'));
cbh(Itime,Jtime).Value = 1;
app.DataLoadingAndFittingTabOutputs.marginalMatrix(Itime,Jtime) = 1;

% set default to marginalize over all non-time data fields
IMarg = Nd+3; %find(contains(rowHeaders,'Marginalize'));
JMarg = find(~contains(colHeaders,'time')&~contains(colHeaders,'Time'));
for i=1:length(JMarg)
    cbh(IMarg,JMarg(i)).Value = 1;
end

% Create Checkboxes
for kR = length(rowHeaders) - 1 : -1 : 1
    for kH = length(colHeaders) : -1 :1
        cbh(kR,kH).Callback = {@marginalizeCheckBoxCallback,[kR,kH],app,cbh,IMarg,Nd};
    end
end

% Create Conditional Fields
for l = 1:length(colHeaders)
    edf(l) = uicontrol(hFig,'Style','edit','Position',[65+140*l 570-(40*(length(rowHeaders)-1)) 120 22],...
        'Callback',{@conditionOnCallback,l,app});
end
% Create Ok Box
applyBox = uicontrol(hFig,'Style','pushbutton','Position',[1100 10 100 40],...
    'BackgroundColor','g','String','Apply','Callback',{@filterAndMarginalize,app},'FontSize',16);
uiwait(hFig)
figure(app.UIFigure);

function conditionOnCallback(hObject,~,fieldId,app)

str = get(hObject,'String');
app.DataLoadingAndFittingTabOutputs.conditionOnArray = string(app.DataLoadingAndFittingTabOutputs.conditionOnArray);
if str
    app.DataLoadingAndFittingTabOutputs.conditionOnArray(fieldId) = str;
else
    app.DataLoadingAndFittingTabOutputs.conditionOnArray(fieldId) = "0";
end

function marginalizeCheckBoxCallback(hObject,~,checkboxId,app,cbh,IMarg,Nd)
value = get(hObject,'Value');
app.DataLoadingAndFittingTabOutputs.marginalMatrix(checkboxId(1),checkboxId(2)) = value;
if checkboxId(1)<=Nd+2
    if value == 0
        cbh(IMarg,checkboxId(2)).Value = 1;
    elseif value == 1
        cbh(IMarg,checkboxId(2)).Value = 0;
    end
elseif checkboxId(1)==IMarg
    if value == 0
    else
        for i=1:Nd+2
            cbh(i,checkboxId(2)).Value = 0;
        end
    end
end
app.DataLoadingAndFittingTabOutputs.marginalMatrix(:) = [cbh.Value];

app.SpeciesForFitPlot.Value = {};
for i=1:Nd 
    if sum(app.DataLoadingAndFittingTabOutputs.marginalMatrix(i,:))
        app.SpeciesForFitPlot.Value(end+1) = app.SpeciesForFitPlot.Items(i);
    end
end

