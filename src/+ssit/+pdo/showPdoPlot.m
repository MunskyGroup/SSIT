function showPdoPlot(app)
% Makes plots of each individual species distortion operator.
arguments
    app
end

ssit.pdo.generatePDO(app);
nSpecies = length(app.FIMTabOutputs.distortionOperator.conditionalPmfs);
figure

% Determine how many plots to make
kPlots = 0;
for i = 1:nSpecies
    Z = app.FIMTabOutputs.distortionOperator.conditionalPmfs{i};
    if min(size(Z))>2
        kPlots=i;
    end
end

% Make plot of PDO
for i = 1:kPlots
    subplot(2,kPlots,i)
    Z = app.FIMTabOutputs.distortionOperator.conditionalPmfs{i};
    if min(size(Z))>2
        contourf(log10(Z))
        xlabel('True Number')
        ylabel('Observed Number')
        title(['Distortion for Species ',num2str(i)])
        set(gca,'fontsize',15)
    end
end
c = colorbar;
set(c.Label,'String','log_{10} p(obs|true)','fontsize',14)

% Make plot of distributions before and after distortion.
if ~isempty(app.FspTabOutputs.solutions)
    px = app.FspTabOutputs.solutions{end}.p;
    py = app.FIMTabOutputs.distortionOperator.computeObservationDist(px);
    Nd = px.dim;
    for i=1:kPlots
        INDS = setdiff([1:Nd],i);
        mdistx = double(px.sumOver(INDS).data);
        mdisty = double(py.sumOver(INDS).data);
        subplot(2,kPlots,kPlots+i) 
        hold off
        plot(mdistx,'linewidth',3)
        hold on
        plot(mdisty,'--','linewidth',3)
        xlabel('Number')
        ylabel('Probability')
        legend({'Actual','Observed'})
        set(gca,'fontsize',15)
    end
end




