%% Make figures of final experiment designs
% fset = TestCases.figureSet;
% f = figure(fset*100+2);
% set(f,'Name',[TestCases.model,' experiment design, ',TestCases.expDesign])
finalExperimentDesign = initialExperiment;
if iModel==8
    for j = 1:min(nRounds-1,5)
        finalExperimentDesign = finalExperimentDesign + TestCases.measurements{j};
    end

else

    for j = 1:nRounds-1
        finalExperimentDesign = finalExperimentDesign + TestCases.measurements{j};
    end
end
% bar(finalExperimentDesign'); hold on;
% set(gca,"FontSize",16)
% xlabel('Measurement Times')
% ylabel('Number of Cells Measured')