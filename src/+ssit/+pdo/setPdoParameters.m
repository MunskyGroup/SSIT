function setPdoParameters(app)
disp('Please use the GUI to enter the distortion operator parameters.')
ssit.parest.PropEditor(app.FIMTabOutputs.FIMPrior,'props',{'Parameter','Value'},1);
app.pdo_parameters_table.Data=[];
