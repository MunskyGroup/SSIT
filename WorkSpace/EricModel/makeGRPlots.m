function makeGRPlots(combinedModel,GRpars)
    combinedGRModel = combinedModel.updateModels(GRpars,false);
    nMods = length(combinedGRModel.SSITModels);
    ModelGroup = cell(nMods,1);
    for i=1:nMods
        %  Update parameters in original models.
        ModelGroup{i} = combinedGRModel.SSITModels{i};
        ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
        ModelGroup{i}.makeFitPlot([],1,[],true,'STD');
    end
end