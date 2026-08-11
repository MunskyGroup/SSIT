function FIM = totalFIM(inputFIMs, observations, covPrior)
    arguments
        inputFIMs
        observations (1, :)
        covPrior = [];
    end

    if isa(observations, "MBExperimentDesign")        
        observations = observations.getAsObservationVector();    
    end

    numConfigurations = size(inputFIMs, 1);
    if numConfigurations ~= length(observations)
        error("totalFIM was called with observations for " + ...
            length(observations) + " configurations but " + ...
            " FIMs for " + numConfigurations + " configurations")
    end

    if isempty(covPrior)
        fimPrior = zeros(size(inputFIMs{1}));
    else
        fimPrior = inv(covPrior);
    end

    numFIMsamples = size(inputFIMs, 2);
    FIM = cell(1,numFIMsamples);

    for sampleIdx = 1:numFIMsamples
        FIM{sampleIdx} = fimPrior;
        for configIdx = 1:numConfigurations
            FIM{sampleIdx} = FIM{sampleIdx} + ...
                observations(configIdx) * inputFIMs{configIdx, sampleIdx};
        end
    end
end