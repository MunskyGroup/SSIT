function [results, Model] = fittingPipelineExample(Model,Args)
% Example pipeline that fits the model to data already contained within the
% model using optional arguments. It runs the pipelne and returns the
% updated model.
arguments
    Model
    Args = struct('maxIter',1000,'display','iter','makePlot',false,'nRounds',1)
end

% Replace missing arguments with defaults.
defaultArgs = struct('maxIter',1000,'display','iter','makePlot',false,'nRounds',1);
argFields = fields(defaultArgs);
for ifield = 1:length(argFields)
    if ~isfield(Args, argFields{ifield})
        Args.(argFields{ifield}) = defaultArgs.(argFields{ifield});
    end
end

% No results returned in this pipeline.
results = [];

% Find initial FSP bounds
[~,~,Model] = Model.solve;  

fitOptions = optimset('Display',Args.display,'MaxIter',Args.maxIter);

% Fit model to data
for iRound = 1:Args.nRounds
    [~,~,~,Model] = Model.maximizeLikelihood([],fitOptions);

    % Find final FSP bounds
    [~,~,Model] = Model.solve;
end

if Args.makePlot
    Model.makeFitPlot;
end

end
