function [results, combinedModel] = multiModelFittingPipelineExample(combinedModel,Args)

% Parse optional arguments
numIter = Args.maxIter;
display = Args.display;
fitOptions = optimset('Display',display,'MaxIter',numIter);

%% Run the fit
Pars = combinedModel.parameters;
Pars = combinedModel.maximizeLikelihood(...
    Pars, fitOptions);
combinedModel = combinedModel.updateModels(Pars,false);

results = [];
