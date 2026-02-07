% This script clears axes in the SSIT GUI
axesNames = {'SsaTrajAxes','SsaHistAxes','FspAxes','SensProbAxes',...
    'SensDerivativeAxes','plotFIMvsTime','FIMEllipseAxes',...
    'PDO_Axis','PDO_Axis2'};

for i=1:length(axesNames)
    cla(app.(axesNames{i}))
end