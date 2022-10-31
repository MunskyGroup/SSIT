function makeMarginalsFsp(app)
% This function creates graphs of the different marginals at the desired
% time point in the FSP 3 Species GUI. This will use the time selected from
% the time slider as the time point from each of the different combinations
% based on the species selected.

T_array = eval(app.FspPrintTimesField.Value);
[~,j] = min(abs(T_array-app.FspTimeSlider.Value));

% Compute the marginal distributions
Nd = app.FspTabOutputs.solutions{j}.p.dim;
if Nd==1
    mdist{1} = double(app.FspTabOutputs.solutions{j}.p.data);
else
    for i=1:Nd
        INDS = setdiff([1:Nd],i);
        mdist{i} = double(app.FspTabOutputs.solutions{j}.p.sumOver(INDS).data);
    end
end

figure()
hold('off')
nSpecies = length(app.ReactionsTabOutputs.varNames);
for iSpecies = 1:nSpecies
    iPlot(iSpecies) = max(contains(app.SpeciestoShowListBoxMargFSP.Value,app.ReactionsTabOutputs.varNames{iSpecies}));
end
%  iPlot = [app.FspMarginalX1CheckBox.Value,app.FspMarginalX2CheckBox.Value,app.FspMarginalX3CheckBox.Value];
cols = {'b', 'r', 'g', 'm', 'c', 'b--', 'r--', 'g--', 'm--', 'c--'};
legs = {};
for i=1:nSpecies
    if iPlot(i)
        stairs([0:length(mdist{i})], [mdist{i};0],cols{i},'LineWidth',2);  hold('on');
        legs{end+1} = char(app.NameTable.Data(i,2));
    end
end
title(sprintf('Marginals at time: t = %1.2f',T_array(j)));
xlabel('Species Count')
ylabel('Probability')
legend(legs);

set(gca,'fontsize',20);
hold('off');

end