function FIMSpecifyPrior(app)           
disp('Please use the GUI to enter the mean and variance for each parameter.')
if isempty(app.FIMTabOutputs.FIMPrior)    
    Npars = size(app.ReactionsTabOutputs.parameters,1);
    app.FIMTabOutputs.FIMPrior = ssit.parest.propsStorage;
    for i = 1:Npars
        fieldName = app.ReactionsTabOutputs.parameters{i,1};
        app.FIMTabOutputs.FIMPrior.props.(fieldName) = ...
            [app.ReactionsTabOutputs.parameters{i,2},0];
    end
end
ssit.parest.PropEditor(app.FIMTabOutputs.FIMPrior,'props',{'Parameter','Mean','Variance'},1);
app.FIMTabOutputs.FIMMatrices = [];