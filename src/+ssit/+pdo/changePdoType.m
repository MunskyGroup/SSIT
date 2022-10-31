function changePdoType(app)
% This function sets the default properties for the user selected
% probabilistic distortion operator.

app.SetDistortionParametersButton.Enable = 1;
app.ShowDistortionPlotButton.Enable = 1;
app.SetDistortionParametersButton_2.Enable = 1;
app.ShowDistortionPlotButton_2.Enable = 1;
app.FIMTabOutputs.distortionOperator = [];
app.FIMTabOutputs.paramsPDO = [];
app.DistortionParameterGuessesandFitsLabel.Visible = 0;
app.pdo_parameters_table.Visible = 0;
app.UpdatePDOParametersButton.Visible = 0;

switch app.DistortionTypeDropDown.Value
    case 'None'
        app.FIMTabOutputs.PDOProperties = [];
        app.SetDistortionParametersButton.Enable = 0;
        app.ShowDistortionPlotButton.Enable = 0;
        app.SetDistortionParametersButton_2.Enable = 0;
        app.ShowDistortionPlotButton_2.Enable = 0;
    case 'Binomial'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS1 = 0.5;
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS2 = 1;
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS3 = 1;
    case 'Binomial - State Dependent'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS1 = '@(x)1 - x/(20+x)';
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS2 = '@(x)1 - x/(20+x)';
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS3 = '@(x)1 - x/(20+x)';
    case 'Binomial - Parametrized'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS1 = '@(x,C)1 - C(1)*x/(20+x)';
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS2 = '@(x,C)1 - C(1)*x/(20+x)';
        app.FIMTabOutputs.PDOProperties.props.CaptureProbabilityS3 = '@(x,C)1 - C(1)*x/(20+x)';
        app.FIMTabOutputs.PDOProperties.props.ParameterGuess = '1';
        app.DistortionParameterGuessesandFitsLabel.Visible = 1;
        app.pdo_parameters_table.Visible = 1;
        app.UpdatePDOParametersButton.Visible = 1;
    case 'Poisson'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS1 = 5;
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS2 = 5;
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS3 = 5;
    case 'Poisson - State Dependent'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS1 = '@(x)0.1*x';
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS2 = '@(x)0.1*x';
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS3 = '@(x)0.1*x';
    case 'Poisson - Parametrized'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS1 = '@(x,C)C(1)*x';
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS2 = '@(x,C)C(1)*x';
        app.FIMTabOutputs.PDOProperties.props.NoiseMeanS3 = '@(x,C)C(1)*x';
        app.FIMTabOutputs.PDOProperties.props.ParameterGuess = '5';
        app.DistortionParameterGuessesandFitsLabel.Visible = 1;
        app.pdo_parameters_table.Visible = 1;
        app.UpdatePDOParametersButton.Visible = 1;
    case 'Binning'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.NumberBinsS1 = 10;
        app.FIMTabOutputs.PDOProperties.props.NumberBinsS2 = 10;
        app.FIMTabOutputs.PDOProperties.props.NumberBinsS3 = 10;
        app.FIMTabOutputs.PDOProperties.props.binTypeS1 = 'lin';
        app.FIMTabOutputs.PDOProperties.props.binTypeS2 = 'lin';
        app.FIMTabOutputs.PDOProperties.props.binTypeS3 = 'lin';
    case 'Custom Function'
        app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
        app.FIMTabOutputs.PDOProperties.props.PDO_S1 = '@(y,x) (y>=x)*5^(y-x)*exp(-5)/factorial(y-x)';
        app.FIMTabOutputs.PDOProperties.props.PDO_S2 = '@(y,x) (y>=x)*5^(y-x)*exp(-5)/factorial(y-x)';
        app.FIMTabOutputs.PDOProperties.props.PDO_S3 = '@(y,x) (y>=x)*5^(y-x)*exp(-5)/factorial(y-x)';
        app.FIMTabOutputs.PDOProperties.props.MaxObservationS1 = 100;
        app.FIMTabOutputs.PDOProperties.props.MaxObservationS2 = 100;
        app.FIMTabOutputs.PDOProperties.props.MaxObservationS3 = 100;
end
