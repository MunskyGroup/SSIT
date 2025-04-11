function [results, combinedModel] = MultiModelFittingPipeline(combinedModel,Args)

results = [];
makePlots = Args.makePlots;
Pars = combinedModel.parameters;

if Args.runFit
    combinedModel = combinedModel.initializeStateSpaces();
    numIter = Args.fitOpts.maxIter;
    display = Args.fitOpts.display;
    fitOptions = optimset('Display',display,'MaxIter',numIter);

    %% Run the fit
    Pars = combinedModel.maximizeLikelihood(Pars, fitOptions);
    combinedModel.parameters = Pars;
    combinedModel = combinedModel.updateModels(Pars,makePlots);
end


%% Run Metropolis Hastings if Requested.
if Args.runMH
    % Compute FIM
    disp('Computing FIM.')
    combinedModel = combinedModel.initializeStateSpaces();
    combinedModel = combinedModel.computeFIMs([],'log');

    disp('Running Metropolis Hastings.')
    MHFitOptions.thin=1;
    MHFitOptions.numberOfSamples=100;
    MHFitOptions.burnIn=0;
    MHFitOptions.progress=true;
    MHFitOptions.numChains = 1;
    MHFitOptions.useFIMforMetHast = true;
    MHFitOptions.saveFile = ['TMPMHRun_',num2str(randi(10000)),'.mat'];
    
    MHfields = fieldnames(Args.MHopts);
    for i = 1:length(MHfields)
        MHFitOptions.(MHfields{i}) = Args.MHopts.(MHfields{i});
    end
    [Pars,~,results.MHResults] = combinedModel.maximizeLikelihood(...
        Pars, MHFitOptions, 'MetropolisHastings');
    
    combinedModel.parameters = Pars;
    combinedModel = combinedModel.updateModels(Pars,makePlots);

    delete(MHFitOptions.saveFile)
end

