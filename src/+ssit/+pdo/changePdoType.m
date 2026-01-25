function changePdoType(app)
% This function sets the default properties for the user selected
% probabilistic distortion operator.

if strcmp(app.DistortionTypeDropDown.Value,'None')
    app.FIMTabOutputs.PDOProperties = [];
    app.SetDistortionParametersButton.Enable = 0;
    app.ShowDistortionPlotButton.Enable = 0;
    app.SpeciesDropDown.Enable = 0;
    app.SolutionTimeDropDown.Enable = 0;
    app.SSITModel.pdoOptions = [];
    cla(app.PDO_Axis);app.PDO_Axis.Visible = false;
    cla(app.PDO_Axis2);app.PDO_Axis2.Visible = false;

else
    speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);
 
    app.FIMTabOutputs.PDOProperties = ssit.parest.propsStorage;
    app.SetDistortionParametersButton.Enable = 1;
    app.ShowDistortionPlotButton.Enable = 1;
    app.FIMTabOutputs.distortionOperator = [];
    app.FIMTabOutputs.paramsPDO = [];
    app.DistortionParameterGuessesandFitsLabel.Visible = 0;
    app.pdo_parameters_table.Visible = 0;
    app.UpdatePDOParametersButton.Visible = 0;
    app.SpeciesDropDown.Enable = 1;
    app.SpeciesDropDown.Items = speciesStochastic;
    app.SolutionTimeDropDown.Enable = 1;
    app.SolutionTimeDropDown.Items = arrayfun(@num2str, app.SSITModel.tSpan, 'UniformOutput', 0);
    app.SolutionTimeDropDown.Value = app.SolutionTimeDropDown.Items(end);
    app.PDO_Axis.Visible = 1;
    app.PDO_Axis2.Visible = 1;


    switch app.DistortionTypeDropDown.Value
        case 'Binomial'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(iSp)]) = 0.5;
            end
        case 'Binomial - State Dependent'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(iSp)]) = '@(x)1 - x/(20+x)';
            end
        case 'Binomial - Parametrized'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['CaptureProbabilityS',num2str(iSp)]) = '@(x,C)1 - C(1)*x/(20+x)';
            end
            app.FIMTabOutputs.PDOProperties.props.ParameterGuess = '1';
            app.DistortionParameterGuessesandFitsLabel.Visible = 1;
            app.pdo_parameters_table.Visible = 1;
            app.UpdatePDOParametersButton.Visible = 1;
        case 'Poisson'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(iSp)]) = 5;
            end
        case 'Poisson - State Dependent'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(iSp)]) = '@(x)0.1*x';
            end
        case 'Poisson - Parametrized'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['NoiseMeanS',num2str(iSp)]) = '@(x,C)C(1)*x';
            end
            app.FIMTabOutputs.PDOProperties.props.ParameterGuess = '5';
            app.DistortionParameterGuessesandFitsLabel.Visible = 1;
            app.pdo_parameters_table.Visible = 1;
            app.UpdatePDOParametersButton.Visible = 1;
        case 'Binning'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['NumberBinsS',num2str(iSp)]) = 10;
                app.FIMTabOutputs.PDOProperties.props.(['binTypeS',num2str(iSp)]) = 'lin';
            end
        case 'Custom Function'
            for iSp = 1:length(speciesStochastic)
                app.FIMTabOutputs.PDOProperties.props.(['PDO_S',num2str(iSp)]) = '@(y,x) (y>=x)*5^(y-x)*exp(-5)/factorial(y-x)';
                app.FIMTabOutputs.PDOProperties.props.(['MaxObservationS',num2str(iSp)]) = 100;
            end
    end
end
updateSpeciesDropBoxes
