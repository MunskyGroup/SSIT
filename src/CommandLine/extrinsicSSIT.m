classdef extrinsicSSIT
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        originalSSIT
        nSamples
        parDistributions
        resultsOriginal
        sampleSSITs
        resultsSSITs
    end

    properties (Dependent)
        averagedResults
    end

    methods
        function obj = extrinsicSSIT(SSIT,parDistributions,nSamples)
            % Construct an instance of this class
            obj.originalSSIT = SSIT;
            obj.nSamples = nSamples;
            obj.parDistributions = parDistributions;
            obj = obj.createAndsolveSampleSSITs();
        end

        function obj = createAndsolveSampleSSITs(obj)
            sampleSSIT = cell(obj.nSamples,1);
            for i = 1:obj.nSamples
                sampleSSIT{i} = obj.originalSSIT;
                sampleSSIT{i}.parameters(:,2) = num2cell(obj.parDistributions());
            end
                  
            obj.originalSSIT.fspOptions.fspTol = 1e-6; 
            [obj.resultsOriginal,bounds] = obj.originalSSIT.solve;
            stateSpace = obj.resultsOriginal.stateSpace;

            resultsStruct = cell(obj.nSamples,1);
            parfor i = 1:obj.nSamples
                sampleSSIT{i}.fspOptions.fspTol = inf;
                sampleSSIT{i}.fspOptions.bounds = bounds;
                resultsStruct{i} = sampleSSIT{i}.solve(stateSpace);
            end

            obj.sampleSSITs = sampleSSIT;
            obj.resultsSSITs = resultsStruct;
        end

        function avResults = get.averagedResults(obj)
            nTime = length(obj.resultsSSITs{1}.fsp);
            avResults.fsp = cell(nTime,1);
            for iTime = 1:nTime
                avResults.fsp{iTime} = obj.resultsSSITs{1}.fsp{iTime};
                avResults.fsp{iTime}.p.data = 0;
                for iSamp = 1:obj.nSamples
                    avResults.fsp{iTime}.p.data = avResults.fsp{iTime}.p.data+...
                       obj.resultsSSITs{iSamp}.fsp{iTime}.p.data/obj.nSamples;
                end
            end
        end

    end
end