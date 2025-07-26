function app = runComputeFim(app,PlotRedundancy)
arguments
    app
    PlotRedundancy=false
end
%RUNCOMPUTEFIM Compute the single-observation Fisher Information Matrix
%across time points.

%% Run Sensitivity calculation for all time points.
% Make a list of all time points, then recalculate Sensitivity matrix if
% necessary (i.e., if there are new times points)
app.FIMTabOutputs.FIMTimes = eval(app.ListofMeasurementTimesEditField.Value);
SensTimes = eval(app.SensPrintTimesEditField.Value);
newTimes = setdiff(app.FIMTabOutputs.FIMTimes,SensTimes);
if ~isempty(newTimes)
    SensTimes = sort([SensTimes,newTimes]);
    app.SensPrintTimesEditField.Value = mat2str(SensTimes);
    app.SensFspTabOutputs.solutions=[];
end
indsFIMTimes = ismember(app.FIMTabOutputs.FIMTimes,SensTimes);

%% Find species that are not observed and will nweed to be summed over in FIM calculation
Nd = size(app.SSITModel.species,1);
indsUnobserved=[];
indsObserved=[];
for i=1:Nd
    if ~contains(app.ObservableSpeciesListBox.Value,app.SSITModel.species{i})
        indsUnobserved=[indsUnobserved,i];
    else
        indsObserved=[indsObserved,i];
    end
end
app.SSITModel.pdoOptions.unobservedSpecies = app.SSITModel.species(indsUnobserved);

%% Run Sensitivity calculation for all parameter samples.
if strcmp(app.ModelUncertaintyDropDown.Value,'MC Sample Over Prior')
    BaseParameters = app.SSITModel.parameters;

    Npars = size(app.ReactionsTabOutputs.parameters,1);
    MN = zeros(Npars,1); VAR = zeros(Npars,1);
    for i = 1:Npars
        fieldName = app.ReactionsTabOutputs.parameters{i,1};
        MN(i) = app.FIMTabOutputs.FIMPrior.props.(fieldName)(1);
        VAR(i) = app.FIMTabOutputs.FIMPrior.props.(fieldName)(2);
    end

    nSamps = app.FIMNumMC.Value+1;
    MHSamples = zeros(nSamps,Npars);
    for k = 1:nSamps
        if k==1
            MHSamples(k,:) = [app.SSITModel.parameters{:,2}];
        else
            MHSamples(k,:) = abs(MN' + sqrt(VAR').*randn(size(MN')));
        end
    end

    [fimResults] = app.SSITModel.computeFIM([],'lin',MHSamples);
    app.FIMTabOutputs.FIMMatrices = fimResults;

    %     if isempty(app.FIMTabOutputs.FIMMatrices)
    %         app.FIMTabOutputs.FIMMatrices = {};
    %         try
    %             f = app.UIFigure;
    %             d_prog_bar = uiprogressdlg(f,'Title','Running MC Sensitivity Computation');
    %             d_prog_bar.Value = 0;
    %         catch
    %             d_prog_bar=[];
    %         end
    %         for k = app.FIMNumMC.Value+1:-1:1
    %             if k==1
    %                 pars = [BaseParameters{:,2}];
    %             else
    %                 pars = MN + sqrt(VAR).*randn(size(MN));
    %             end
    %             for j = 1:Npars
    %                 app.ReactionsTabOutputs.parameters{j,2} = pars(j);
    %             end
    %             app = runSensitivity(app,false);
    %             sensoutputs{k} = app.SensFspTabOutputs.solutions.data(indsFIMTimes);
    %
    %             ssit.pdo.generatePDO(app,[],sensoutputs{k},indsObserved);
    %
    %             for it=length(sensoutputs{k}):-1:1
    % %                 try
    %                     if isempty(indsUnobserved)
    %                         F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p, sensoutputs{it}.S, app.FIMTabOutputs.distortionOperator);
    %                     else
    %                         % Remove unobservable species.
    %                         redS = sensoutputs{it}.S;
    %                         for ir = 1:length(redS)
    %                             redS(ir) = sensoutputs{it}.S(ir).sumOver(indsUnobserved);
    %                         end
    %                         F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p.sumOver(indsUnobserved), redS, app.FIMTabOutputs.distortionOperator);
    %                     end
    % %                 catch
    % %                     ssit.pdo.generatePDO(app,[],sensoutputs{k});
    % %                     app.FIMTabOutputs.distortionOperator.conditionalPmfs =  app.FIMTabOutputs.distortionOperator.conditionalPmfs(indsObserved);
    % %                     if isempty(indsUnobserved)
    % %                         F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p, sensoutputs{it}.S, app.FIMTabOutputs.distortionOperator);
    % %                     else
    % %                         % Remove unobservable species.
    % %                         redS = sensoutputs{it}.S;
    % %                         for ir = 1:length(redS)
    % %                             redS(ir) = redS(ir).sumOver(indsUnobserved);
    % %                         end
    % %                         F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p.sumOver(indsUnobserved), redS, app.FIMTabOutputs.distortionOperator);
    % %                     end
    % %                 end
    %                 app.FIMTabOutputs.FIMMatrices{it,k} = F;
    %             end
    %             d_prog_bar.Value = (app.FIMNumMC.Value+1-k+1)/(app.FIMNumMC.Value+1);
    %         end
    %     end
else

    if (isempty(app.SensFspTabOutputs.solutions))
        app = runSensitivity(app);
    end
    sensoutputs = app.SensFspTabOutputs.solutions;

    % Call function to generate the distortion operator.
    % ssit.pdo.generatePDO(app,[],sensoutputs,indsObserved);


    %% Compute FIM for every time point.
    [fimResults] = app.SSITModel.computeFIM(sensoutputs,'lin');
    app.FIMTabOutputs.FIMMatrices = fimResults;
    %     for it=length(sensoutputs):-1:1
    % %         try
    %             if isempty(indsUnobserved)
    %                 F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p, sensoutputs{it}.S, app.FIMTabOutputs.distortionOperator);
    %             else
    %                 % Remove unobservable species.
    %                 redS = sensoutputs{it}.S;
    %                 for ir = 1:length(redS)
    %                     redS(ir) = sensoutputs{it}.S(ir).sumOver(indsUnobserved);
    %                 end
    %                 F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p.sumOver(indsUnobserved), redS, app.FIMTabOutputs.distortionOperator);
    %             end
    % %         catch
    % %             ssit.pdo.generatePDO(app,[],sensoutputs);
    % %             app.FIMTabOutputs.distortionOperator.conditionalPmfs =  app.FIMTabOutputs.distortionOperator.conditionalPmfs(indsObserved);
    % %             if isempty(indsUnobserved)
    % %                 F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p, sensoutputs{it}.S, app.FIMTabOutputs.distortionOperator);
    % %             else
    % %                 % Remove unobservable species.
    % %                 redS = sensoutputs{it}.S;
    % %                 for ir = 1:length(redS)
    % %                     redS(ir) = redS(ir).sumOver(indsUnobserved);
    % %                 end
    % %                 F = ssit.fim.computeSingleCellFim(sensoutputs{it}.p.sumOver(indsUnobserved), redS, app.FIMTabOutputs.distortionOperator);
    % %             end
    % %         end
    %         app.FIMTabOutputs.FIMMatrices{it,1} = F;
% end
end

switch app.FIMMetricorParameterDropDown.Value
    case 'Determinant'
        met = @(A)det(A);
        app.plotFIMvsTime.YLabel.String = 'det(FIM)';
    case 'Smallest Eigenvalue'
        met = @(A)min(eig(A));
        app.plotFIMvsTime.YLabel.String = 'min(\lambda_{FIM})';
    case 'Trace'
        met = @(A)trace(A);
        app.plotFIMvsTime.YLabel.String = 'trace(FIM)';
    otherwise
        k = find(strcmp(app.SensParDropDown.Items,app.FIMMetricorParameterDropDown.Value));
        ek = zeros(1,length(app.FIMTabOutputs.FIMMatrices{1,1}));ek(k) = 1;
        met = @(A)(-ek*inv(A)*ek');
        app.plotFIMvsTime.YLabel.String = ['variance reduction for ',app.FIMMetricorParameterDropDown.Value];
end
NSamp = size(app.FIMTabOutputs.FIMMatrices,2);
NT = size(app.FIMTabOutputs.FIMMatrices,1);
FIMMetricVsTime = zeros(NT,NSamp);
for k = 1:NSamp
    FIMTot = 0*app.FIMTabOutputs.FIMMatrices{1,1};
    for it=1:NT
        FIMTot = FIMTot+app.FIMTabOutputs.FIMMatrices{it,k};
    end

    for it=1:NT
        A = app.FIMTabOutputs.FIMMatrices{it,k};
        if PlotRedundancy
            FIMMetricVsTime(it,k) = met(FIMTot) - met(FIMTot-A);
        else
            FIMMetricVsTime(it,k) = met(A);
        end
    end
end

if rank(FIMTot)<size(FIMTot,1)
    warndlg(['FIM has rank r=',num2str(rank(FIMTot)),' - You cannot identify ',num2str(size(FIMTot,1)),' parameters with this experiment.'],'Rank Warning')
    cla(app.plotFIMvsTime);
    return
end


if NSamp==1
    plot(app.plotFIMvsTime,app.FIMTabOutputs.FIMTimes,FIMMetricVsTime,'o','linewidth',3)
    plot(app.FIMTabOutputs.FIMTimes,FIMMetricVsTime,'o','linewidth',3)
else
    errorbar(app.plotFIMvsTime,app.FIMTabOutputs.FIMTimes,mean(FIMMetricVsTime,2),std(FIMMetricVsTime,0,2),'o','linewidth',3)
    errorbar(app.FIMTabOutputs.FIMTimes,mean(FIMMetricVsTime,2),std(FIMMetricVsTime,0,2),'o','linewidth',3)
end
app.plotFIMvsTime.YLim = [min(FIMMetricVsTime(:)),max(FIMMetricVsTime(:))];

set(gca,'ylim',[min(FIMMetricVsTime(:)),max(FIMMetricVsTime(:))],...
    'fontsize',16);
xlabel('time')
switch app.FIMMetricorParameterDropDown.Value
    case 'Determinant'
        ylabel('det(FIM)');
    case 'Smallest Eigenvalue'
        ylabel('min(\lambda_{FIM})');
    case 'Trace'
        ylabel('trace(FIM)');
    otherwise
        ylabel(['variance reduction for ',app.FIMMetricorParameterDropDown.Value]);
end

app.FIMParameter1.Items = app.SSITModel.parameters(:,1);
app.FIMParameter2.Items = app.SSITModel.parameters(:,1);
app.FIMParameter2.Value = app.SSITModel.parameters(2,1);



