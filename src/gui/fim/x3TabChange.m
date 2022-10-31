function x3TabChange(app)
%% This Function prompts the user to pick type of noise error in image.
% This prompt is in the FIM tab and located above the Uniform Processing 
% Button in a drop down menu. Four options will be given. Error Free, 
% Binomial, Poisson, and Gaussian.
value = app.ImageProcessingErrorsDropDown_5.Value;
if strcmp(value,'Error Free')
    app.ProbEditField_3.Enable = 'off';
    app.MeanEditField_3.Enable = 'off';
    app.VarianceEditField_3.Enable = 'off';
    app.ProbiLabel_3.Text = 'Prob. @(i)';
    app.ProbEditField_3.Visible = 'on';
    app.MeanEditField_3.Visible = 'off';
elseif strcmp(value,'Binomial')
    app.ProbEditField_3.Enable = 'on';
    app.MeanEditField_3.Enable = 'off';
    app.VarianceEditField_3.Enable = 'off';
    app.MeanEditField_3.Visible = 'off';
    app.ProbEditField_3.Visible = 'on';
    app.ProbiLabel_3.Text = 'Prob. @(i)';
elseif strcmp(value,'Poisson')
    app.ProbEditField_3.Enable = 'off';
    app.MeanEditField_3.Enable = 'on';
    app.VarianceEditField_3.Enable = 'off';
    app.MeanEditField_3.Visible = 'on';
    app.ProbEditField_3.Visible = 'off';
    app.ProbiLabel_3.Text = 'Mean @(i)';
elseif strcmp(value,'Gaussian')
    app.ProbEditField_3.Enable = 'off';
    app.MeanEditField_3.Enable = 'on';
    app.VarianceEditField_3.Enable = 'on';
    app.MeanEditField_3.Visible = 'on';
    app.ProbEditField_3.Visible = 'off';
    app.ProbiLabel_3.Text = 'Mean @(i)';
end
    