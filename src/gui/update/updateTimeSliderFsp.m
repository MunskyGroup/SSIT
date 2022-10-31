function updateTimeSliderFsp(app,event)
% Updates the time slider when the time array is changed to match on the graphs
arguments
    app
    event = [];
end
if ~isempty(event)
    valueStr = event.Source.Value;
else
    valueStr = app.PrintTimesEditField.Value;
end
value = eval(valueStr);

app.PrintTimesEditField.Value = valueStr;
app.FspPrintTimesField.Value = valueStr;
app.SensPrintTimesEditField.Value = valueStr;
app.ListofMeasurementTimesEditField.Value = valueStr;

app.FspTimeSlider.Limits = [min(value),max(value)];
app.FspTimeSlider.Value = max(value);
L = length(value);
if L<4
    app.FspTimeSlider.MajorTicks = value;
else
    app.FspTimeSlider.MajorTicks = value(1:floor(L/4):L);
end

app.SensPlotTimeSlider.Limits = [min(value),max(value)];
app.SensPlotTimeSlider.Value = max(value);
L = length(value);
if L<4
    app.SensPlotTimeSlider.MajorTicks = value;
else
    app.SensPlotTimeSlider.MajorTicks = value(1:floor(L/4):L);
end

app.SsaTimeSlider.Limits = [min(value),max(value)];
L = length(value);
if L<4
    app.SsaTimeSlider.MajorTicks = value;
else
    app.SsaTimeSlider.MajorTicks = value(1:floor(L/4):L);
end

app.FIMTabOutputs.CellsPerTimePoint=[];
app.FIMTabOutputs.FIMMatrices=[];
app.FIMTabOutputs.NcOptimized=[];
app.SensFspTabOutputs.solutions=[];

end
