addpath(genpath('../src'));
%% Call Pipeline to Fit Model

% Specify pipeline to apply to model and arguments
% ("../../SSIT/src/exampleData/examplePipelines/fittingPipelineExample.m") 
Pipeline = 'fittingPipelineExample';
pipelineArgs.maxIter = 1000;
pipelineArgs.display = 'iter';
pipelineArgs.makePlot = false;
pipelineArgs.nRounds = 1;

%% Launch cluster jobs for all genes.
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
TAB = readtable(DataFileName);
geneNames = fields(TAB);

for iRound = 1:5
    for iGene = 2:length(geneNames)-4
        modelName = ['Model_',geneNames{iGene}];
        saveName = ['seqModels/',modelName];
        logfile = ['logFiles/log',modelName];
        if iGene==2
            load(saveName)
            eval([modelName,'.formPropensitiesGeneral;']);
        end
        cmd = SSIT.generateCommandLinePipeline(saveName,modelName,Pipeline,...
            pipelineArgs=pipelineArgs,saveFileOut=saveName, ...
            logFile=logfile,runNow=true,runOnCluster=true);
        pause(60); % Pause to allow for matlab licenses to reset.
    end
end
