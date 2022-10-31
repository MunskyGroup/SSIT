function UniformProcessingX2(app)
%% This function computes uniform processing for the x2 tab.
% This functions activates when the Uniform Processing button is pushed.
ImageProcessValue = app.ImageProcessingErrorsDropDown_4.Value;
ProbabilityValue = app.ProbEditField_2.Value;
MeanValue = app.MeanEditField_2.Value;
VarianceValue = app.VarianceEditField_2.Value;
app.ImageProcessingErrorsDropDown.Value = ImageProcessValue;
app.ProbEditField.Value = ProbabilityValue;
app.MeanEditField.Value = MeanValue;
app.VarianceEditField.Value = VarianceValue;
app.ImageProcessingErrorsDropDown_5.Value = ImageProcessValue;
app.ProbEditField_3.Value = ProbabilityValue;
app.MeanEditField_3.Value = MeanValue;
app.VarianceEditField_3.Value = VarianceValue;
x1TabChange(app)
x3TabChange(app)
end