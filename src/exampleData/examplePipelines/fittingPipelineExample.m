function [results, Model] = fittingPipelineExample(Model,Args)
results = [];

% Parse optional arguments
numIter = Args.maxIter;
display = Args.display;
makePlot = Args.makePlot;

% Find initial FSP bounds
[~,~,Model] = Model.solve;  

fitOptions = optimset('Display',display,'MaxIter',numIter);

% Fit model to data
[~,~,~,Model] = Model.maximizeLikelihood([],fitOptions);

% Find final FSP bounds
[~,~,Model] = Model.solve;  

if makePlot
    Model.makeFitPlot;
end

end
