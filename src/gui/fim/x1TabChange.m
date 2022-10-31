function x1TabChange(app)
%% This Function prompts the user to pick type of noise error in image.
% This prompt is in the FIM tab and located above the Uniform Processing 
% Button in a drop down menu. Four options will be given. Error Free, 
% Binomial, Poisson, and Gaussian.
value = app.ImageProcessingErrorsDropDown.Value;
if strcmp(value,'Error Free')
    app.ProbEditField.Enable = 'off';
    app.MeanEditField.Enable = 'off';
    app.VarianceEditField.Enable = 'off';
    app.ProbiLabel.Text = 'Prob. @(i)';
    app.ProbEditField.Visible = 'on';
    app.MeanEditField.Visible = 'off';
elseif strcmp(value,'Binomial')
    app.ProbEditField.Enable = 'on';
    app.MeanEditField.Enable = 'off';
    app.VarianceEditField.Enable = 'off';
    app.MeanEditField.Visible = 'off';
    app.ProbEditField.Visible = 'on';
    app.ProbiLabel.Text = 'Prob. @(i)';
elseif strcmp(value,'Poisson')
    app.ProbEditField.Enable = 'off';
    app.MeanEditField.Enable = 'on';
    app.VarianceEditField.Enable = 'off';
    app.MeanEditField.Visible = 'on';
    app.ProbEditField.Visible = 'off';
    app.ProbiLabel.Text = 'Mean @(i)';
elseif strcmp(value,'Gaussian')
    app.ProbEditField.Enable = 'off';
    app.MeanEditField.Enable = 'on';
    app.VarianceEditField.Enable = 'on';
    app.MeanEditField.Visible = 'on';
    app.ProbEditField.Visible = 'off';
    app.ProbiLabel.Text = 'Mean @(i)';
end
    