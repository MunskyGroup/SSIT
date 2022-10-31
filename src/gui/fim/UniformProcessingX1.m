function UniformProcessingX1(app)
%% This function computes uniform processing for the x1 tab.
% This functions activates when the Uniform Processing button is pushed.
ImageProcessValue = app.ImageProcessingErrorsDropDown.Value;
ProbabilityValue = app.ProbEditField.Value;
MeanValue = app.MeanEditField.Value;
VarianceValue = app.VarianceEditField.Value;
app.ImageProcessingErrorsDropDown_4.Value = ImageProcessValue;
app.ProbEditField_2.Value = ProbabilityValue;
app.MeanEditField_2.Value = MeanValue;
app.VarianceEditField_2.Value = VarianceValue;
app.ImageProcessingErrorsDropDown_5.Value = ImageProcessValue;
app.ProbEditField_3.Value = ProbabilityValue;
app.MeanEditField_3.Value = MeanValue;
app.VarianceEditField_3.Value = VarianceValue;
x2TabChange(app)
x3TabChange(app)
end