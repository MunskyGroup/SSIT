function FIMEstimateUncertainty(app)

parNames = app.SSITModel.parameters(:,1);
parIndices = containers.Map(parNames, 1:length(parNames));

iX = parIndices(app.FIMParameter1.Value);
iY = parIndices(app.FIMParameter2.Value);

NSamp = size(app.FIMTabOutputs.FIMMatrices,2);
NT = size(app.FIMTabOutputs.FIMMatrices,1);

% figure
hold(app.FIMEllipseAxes,'off')
for js = 1:NSamp
    %% Compute Total FIM
    FIM = zeros(length(parIndices));
    for it = 1:NT
        fieldName = num2str(it);
        Nc = app.FIMTabOutputs.CellsPerTimePoint.props.(['t',fieldName])(2);
        fim = app.FIMTabOutputs.FIMMatrices{it,js};
        FIM = FIM + Nc*fim;
    end
    % covMatrix = inv(FIM([iX iY], [iX iY]));
    covMat = inv(FIM);
    covMatrix = covMat([iX iY], [iX iY]);
    mu = [app.ReactionsTabOutputs.parameters{iX, 2}, app.ReactionsTabOutputs.parameters{iY, 2}]';
    if app.FIMLogTransformCheckBox.Value
        covMatrix = ((1./mu)*(1./mu)').*covMatrix;
        mu = log(mu);
    end

    if iX==iY
        sig = sqrt(covMatrix(1,1));
        tt = linspace(mu(1) - 4*sig,mu(1) + 4*sig);
        pdf = normpdf(tt,mu(1),sig);
        plot(app.FIMEllipseAxes,tt,pdf);
        hold(app.FIMEllipseAxes,'on')
    else
        [ellXs, ellYs] = ssit.parest.ellipse(mu,covMatrix);
        plot(app.FIMEllipseAxes, ellXs, ellYs,'linewidth',3);
        hold(app.FIMEllipseAxes,'on')
        % plot(ellXs, ellYs,'linewidth',3);
        hold('on')
    end
end

set(gca,'fontsize',16)
if iX==iY
    app.FIMEllipseAxes.Title.String=['Expected Uncertainty in ',app.FIMParameter1.Value];
    if app.FIMLogTransformCheckBox.Value
        app.FIMEllipseAxes.XLabel.String = ['log_{10} (',app.FIMParameter1.Value,')'];
    else
        app.FIMEllipseAxes.XLabel.String = [app.FIMParameter1.Value];
        app.FIMEllipseAxes.YLabel.String = 'Probability';
    end
    ylabel('Probability')
else
    app.FIMEllipseAxes.Title.String=['Expected Joint Uncertainty in ',app.FIMParameter1.Value,' and ',app.FIMParameter2.Value];
    if app.FIMLogTransformCheckBox.Value
        app.FIMEllipseAxes.XLabel.String = ['log_{10} (',app.FIMParameter1.Value,')'];
        app.FIMEllipseAxes.YLabel.String = ['log_{10} (',app.FIMParameter2.Value,')'];
    else
        app.FIMEllipseAxes.XLabel.String = ['Parameter: ',app.FIMParameter1.Value];
        app.FIMEllipseAxes.YLabel.String = ['Parameter: ',app.FIMParameter2.Value];
    end
end

end