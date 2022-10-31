function FIMManuallyAllocateMeasurements(app)
if isempty(app.FIMTabOutputs.CellsPerTimePoint)    
    app.FIMTabOutputs.CellsPerTimePoint = ssit.parest.propsStorage;
    for i = 1:length(app.FIMTabOutputs.FIMTimes)
        fieldName = num2str(i);
        app.FIMTabOutputs.CellsPerTimePoint.props.(['t',fieldName]) = ...
            [app.FIMTabOutputs.FIMTimes(i),100];
    end
end
ssit.parest.PropEditor(app.FIMTabOutputs.CellsPerTimePoint,'props',{'_','time (do not edit)','Number of Cells'},1);