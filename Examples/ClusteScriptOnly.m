addpath(genpath('../src'));
%% Call Pipeline to Fit Model

% Specify pipeline to apply to model and arguments
% ("../../SSIT/src/exampleData/examplePipelines/fittingPipelineExample.m") 
Pipeline = 'fittingPipelineExample';
pipelineArgs.maxIter = 1000;
pipelineArgs.display = 'iter';
pipelineArgs.makePlot = false;
pipelineArgs.nRounds = 5;

%% Launch cluster jobs for all genes.
DataFileName = 'data/Raw_DEX_UpRegulatedGenes_ForSSIT.csv';
TAB = readtable(DataFileName);
geneNames = fields(TAB);

for iGene = 1:length(geneNames)-4
    modelName = ['Model_',geneNames{iGene}];
    saveName = ['seqModels/',modelName];    
    logfile = ['logFiles/log',modelName];
    cmd = SSIT.generateCommandLinePipeline(saveName,modelName,[],Pipeline,...
        pipelineArgs,saveName,logfile,1,1);
    pause(0.2);
end