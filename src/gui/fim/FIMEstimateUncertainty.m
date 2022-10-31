function FIMEstimateUncertainty(app)

parNames = app.Model.parameterNames;
parIndices = containers.Map(parNames, 1:length(parNames));

iX = parIndices(app.FIMSpecies1.Value);
iY = parIndices(app.FIMSpecies2.Value);

NSamp = size(app.FIMTabOutputs.FIMMatrices,2);
NT = size(app.FIMTabOutputs.FIMMatrices,1);

figure
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
    
    mu = [app.ReactionsTabOutputs.parameters{iX, 2}, app.ReactionsTabOutputs.parameters{iY, 2}]';
    covMatrix = inv(FIM([iX iY], [iX iY]));
    
    if app.FIMLogTransformCheckBox.Value
        covMatrix = ((1./mu)*(1./mu)').*covMatrix;
        mu = log(mu);
    end
    
    [ellXs, ellYs] = ssit.parest.ellipse(mu,covMatrix);
    plot(app.FIMEllipseAxes, ellXs, ellYs,'linewidth',3);
    hold(app.FIMEllipseAxes,'on')
    plot(ellXs, ellYs,'linewidth',3);
    hold('on')
end
set(gca,'fontsize',16)
title(['Expected Joint Uncertainty in ',app.FIMSpecies1.Value,' and ',app.FIMSpecies2.Value])
xlabel(['Parameter: ',app.FIMSpecies1.Value])
ylabel(['Parameter: ',app.FIMSpecies2.Value])

end