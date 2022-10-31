function FIMTableEdit(app)
%% This code will prompt the users to input a time frame.
% The input must be in [start:increment:end]
% For example [0:0.1:100]

value = app.TimePointsEditField.Value;
TimePointData = str2num(value);
CellNumData = app.UITable6.Data(:,2);
app.UITable6.Data = [];
app.UITable6.Data(:,1) = TimePointData;

if numel(TimePointData) < numel(CellNumData)
   app.UITable6.Data(:,2) = CellNumData(1:numel(TimePointData));
elseif numel(TimePointData) > numel(CellNumData)
    app.UITable6.Data(1:numel(CellNumData),2) = CellNumData;
    app.UITable6.Data(numel(CellNumData)+1:numel(TimePointData),2) = 100*...
    ones(numel(TimePointData)-numel(CellNumData),1);
end
end