function showPdoPlot(app)
% Makes plots of each individual species distortion operator.
arguments
    app
end

if isempty(app.SSITModel.pdoOptions.PDO)
    [~,app.SSITModel.pdoOptions.PDO] = ssit.pdo.generatePDO(app);
end
% nSpecies = length(app.FIMTabOutputs.distortionOperator.conditionalPmfs);
% figure

% Determine how many plots to make
% kPlots = 0;
% for i = 1:nSpecies
iSp = find(strcmp(app.SSITModel.species,app.SpeciesDropDown.Value));
Z = app.SSITModel.pdoOptions.PDO.conditionalPmfs{iSp};
% if min(size(Z))>2
% kPlots=i;
% end
% end

% Make plot of PDO
% for i = 1:kPlots
% subplot(2,kPlots,i)
% Z = app.FIMTabOutputs.distortionOperator.conditionalPmfs{i};
% if min(size(Z))>2
contourf(app.PDO_Axis,log10(Z))
% xlabel('True Number')
% ylabel('Observed Number')
% title(['Distortion for Species ',num2str(i)])
% set(gca,'fontsize',15)
% end
% end
% c = colorbar;
% set(c.Label,'String','log_{10} p(obs|true)','fontsize',14)

% Make plot of distributions before and after distortion.
% if ~isempty(app.FspTabOutputs.solutions)
iTime = find(app.SSITModel.tSpan==eval(app.SolutionTimeDropDown.Value));
px = app.FspTabOutputs.solutions{iTime}.p;
py = app.SSITModel.pdoOptions.PDO.computeObservationDist(px);
Nd = px.dim;
% for i=1:kPlots
INDS = setdiff([1:Nd],iSp);
if isempty(INDS)
    mdistx = double(px.data);
    mdisty = double(py.data);
else
    mdistx = double(px.sumOver(INDS).data);
    mdisty = double(py.sumOver(INDS).data);
end
% subplot(2,kPlots,kPlots+i)
% hold off
hold(app.PDO_Axis2,"off")
plot(app.PDO_Axis2,mdistx,'linewidth',3)
hold(app.PDO_Axis2,"on")
plot(app.PDO_Axis2,mdisty,'--','linewidth',3)
% xlabel('Number')
% ylabel('Probability')
legend(app.PDO_Axis2,{'True','Observed'})
% set(gca,'fontsize',15)
% end
% end




