function UniformProcessingX3(app)
%% This function computes uniform processing for the x3 tab.
% This functions activates when the Uniform Processing button is pushed.
ImageProcessValue = app.ImageProcessingErrorsDropDown_5.Value;
ProbabilityValue = app.ProbEditField_3.Value;
MeanValue = app.MeanEditField_3.Value;
VarianceValue = app.VarianceEditField_3.Value;
app.ImageProcessingErrorsDropDown_4.Value = ImageProcessValue;
app.ProbEditField_2.Value = ProbabilityValue;
app.MeanEditField_2.Value = MeanValue;
app.VarianceEditField_2.Value = VarianceValue;
app.ImageProcessingErrorsDropDown.Value = ImageProcessValue;
app.ProbEditField.Value = ProbabilityValue;
app.MeanEditField.Value = MeanValue;
app.VarianceEditField.Value = VarianceValue;
x1TabChange(app)
x2TabChange(app)
end