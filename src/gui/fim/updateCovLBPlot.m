function app = updateCovLBPlot(app)
%UPDATEELLIPSEPLOT Summary of this function goes here
%   Detailed explanation goes here

% applyLogTransform = app.FIMLogTransformCheckBox.Value;
parNames = app.Model.parameterNames;
parIndices = containers.Map(parNames, 1:length(parNames));

% timeIdx = ceil(app.FIMTimeIndexEditField.Value);
iX = parIndices(app.FIMSpecies1.Value);
iY = parIndices(app.FIMSpecies2.Value);

mu = [app.ReactionsTabOutputs.parameters{iX, 2}, app.ReactionsTabOutputs.parameters{iY, 2}]';
covMatrix = inv(app.FIMTabOutputs.fimsSingle{timeIdx}([iX iY], [iX iY]));

% if applyLogTransform
%     covMatrix = ((1./mu)*(1./mu)').*covMatrix;
%     mu = log(mu);    
% end

[ellXs, ellYs] = ssit.parest.ellipse(mu,covMatrix);

plot(app.FIMEllipseAxes, ellXs, ellYs);

end

