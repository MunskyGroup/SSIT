function [results, Model] = fittingPipelineExample(Model,Args)

% Parse optional arguments
numIter = Args.maxIter;
display = Args.display;

% Find initial FSP bounds
[~,Model.fspOptions.bounds] = Model.solve;  

fitOptions = optimset('Display',display,'MaxIter',numIter);

% Fit model to data
results.Modelpars = [Model.parameters{:,2}];
[results.Modelpars,results.Model_likelihood] = Model.maximizeLikelihood(results.Modelpars,...
                                                            fitOptions);

% Update Model Parameters
Model.parameters(:,2) = num2cell(results.Modelpars);

% Find final FSP bounds
[~,Model.fspOptions.bounds] = Model.solve;  

end
