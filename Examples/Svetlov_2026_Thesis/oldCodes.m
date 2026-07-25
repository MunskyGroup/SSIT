% make figures if requested.
if showPlots
    % Plotting options
    mhPlotScale = 'log10';  % Show MH and FIM plots in log10 scale.

    FIMOptNextExptReduced = cell(size(FIMOptNextExpt));
    for i = 1:nFIMsamples
        FIMOptNextExptReduced{i} = FIMOptNextExpt{i}(ModelGuess{1}.fittingOptions.modelVarsToFit,...
            ModelGuess{1}.fittingOptions.modelVarsToFit); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    f = figure;
    f.Name = ['Current MH Results and Next FIM Prediction (Round ',num2str(iExpt),')'];
    ModelGuess{1}.plotMHResults(MHResults,FIMOptNextExptReduced,fimScale,mhPlotScale,f)

    if iExpt>1
        f = figure;
        f.Name = ['Current MH Results and Previous FIM Prediction (Round ',num2str(iExpt),')'];
        ModelGuess{1}.plotMHResults(MHResults,FIMpredNextExpt{iExpt-1},fimScale,mhPlotScale,f)
    end

    f = figure;
    f.Name = ['Current MH Results and Perfect FIM Prediction (Round ',num2str(iExpt),')'];
    ModelTrue{1}.plotMHResults(MHResults,FIMCurrentExpt_True,fimScale,mhPlotScale,f)

    figure;
    title(['Number of cells Measured at each Time Point (Round ',num2str(iExpt),')']);
    bar(ModelGuess{1}.tSpan,nTotalCells,0.4, 'stacked')
    ylabel('Number of Cells Measured');
    xlabel('time [min]');
    ylim([0 300]);
end

% Save results
save(saveExpName,'parametersFound','FIMcurrentExptSaved','FIMcurrentExptTrueSaved','covMH',...
    'covLogMH','exptDesigns','MHResultsSaved','FIMpredNextExpt')