function showPdoPlot(app)
% Makes plots of each individual species distortion operator.
arguments
    app
end

if isempty(app.SSITModel.pdoOptions.PDO)
    [~,app.SSITModel.pdoOptions.PDO] = ssit.pdo.generatePDO(app);
end

% Determine how many plots to make
indSp = find(strcmp(app.SpeciesDropDown.Items,app.SpeciesDropDown.Value));
Z = app.SSITModel.pdoOptions.PDO.conditionalPmfs{indSp};
app.PDO_Axis.Visible = true;
contourf(app.PDO_Axis,log10(Z))

% Make plot of distributions before and after distortion.
% if ~isempty(app.FspTabOutputs.solutions)
iTime = find(app.SSITModel.tSpan==eval(app.SolutionTimeDropDown.Value));
if ~isempty(app.SSITModel.Solutions)||~isfield(app.SSITModel.Solutions,'fsp')
    app.SSITModel.solutionScheme='fsp';
    [~,~,app.SSITModel] = app.SSITModel.solve;
end
%%
maxNum = app.SSITModel.Solutions.fsp{iTime}.p.data.size;
kSp = 0;

speciesStochastic = setdiff(app.SSITModel.species,app.SSITModel.hybridOptions.upstreamODEs);


for iS = 1:length(speciesStochastic)
    if max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic(iS)))
        maxNum(iS) = 0;
        curNum(iS) = 0;
    else
        kSp = kSp+1;
        curNum(iS) = size(app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp},2);
    end
end
if max(maxNum-curNum)>0
    [~,app.SSITModel] = app.SSITModel.generatePDO([],[],[],[],maxNum);
end

Nd = app.SSITModel.Solutions.fsp{iTime}.p.dim;
kSp = 0;
for iSp = 1:length(speciesStochastic)
    if isempty(app.SSITModel.pdoOptions.unobservedSpecies)||~max(strcmpi(app.SSITModel.pdoOptions.unobservedSpecies,speciesStochastic{iSp}))
        kSp = kSp+1;
        if kSp==indSp
            INDS = setdiff([1:Nd],iSp);
            if ~isempty(INDS)
                px = double(app.SSITModel.Solutions.fsp{iTime}.p.sumOver(INDS).data);
            else
                px = double(app.SSITModel.Solutions.fsp{iTime}.p.data);
            end
            py = app.SSITModel.pdoOptions.PDO.conditionalPmfs{kSp}*px;
            break
        end
    end
end




%%
app.PDO_Axis2.Visible = true;
hold(app.PDO_Axis2,"off")
plot(app.PDO_Axis2,px,'linewidth',3)
hold(app.PDO_Axis2,"on")
plot(app.PDO_Axis2,py,'--','linewidth',3)
% xlabel('Number')
% ylabel('Probability')
legend(app.PDO_Axis2,{'True','Observed'})
% set(gca,'fontsize',15)
% end
% end




