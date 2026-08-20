function [logL] = getLogLikelihoodOfDataGivenModels(parametersLogSpace, ...
    models, modelsHaveData, stateSpaces)
%GETLOGLIKELIHOODOFDATAGIVENMODELS This function takes a group of models
% that share a parameter set and then calculates the log-likelihood of 
% available data given those models, taking into account (if it exists) a 
% common prior.
arguments (Input)
    parametersLogSpace
    models (1, :) %SSIT Debugging: this is a cell in the CDC code
    modelsHaveData (1, :) logical
    stateSpaces (1, :)
end

arguments (Output)
    logL (1, 1) double {mustBeNonpositive}   
end

if isempty(stateSpaces)
    stateSpaces = cell(1, length(models));
end

logL = 0;

% We only want to include the prior in one likelihood calculation;
% otherwise, it will be amplified in subsequent data sets.

modelIdxForIncludingPrior = 0;
if any(modelsHaveData)
    % If at least one model has data, sort the list in ascending order; the
    % last index will be guaranteed to correspond to a "true," i.e. a model
    % with data.

    [~, indices] = sort(modelsHaveData);
    modelIdxForIncludingPrior = modelsHaveData(indices(end));
end

% DEBUGGING: Unrecognized function or variable 'err_loc' in some workers
% when parallelizing this using a hybrid CDC/MB code.

for modelIdx = 1:length(models)
    % If there are no data associated with the model, there is no
    % corresponding likelihood to calculate:

    if modelsHaveData(modelIdx)
        if modelIdx == modelIdxForIncludingPrior
            % DEBUGGING:
            logL = logL + models{modelIdx}.computeLikelihood(...
                exp(parametersLogSpace), stateSpaces{modelIdx}); 
            % logL = logL + models(modelIdx).computeLikelihood(...
            %     exp(parametersLogSpace), stateSpaces{modelIdx});        
        else
            % DEBUGGING:
            tempModel = models{modelIdx}; 
            % tempModel = models(modelIdx); 
            tempModel.fittingOptions.logPrior = [];
            logL = logL + tempModel.computeLikelihood(...
                exp(parametersLogSpace), stateSpaces{modelIdx});            
        end
    end
end

end