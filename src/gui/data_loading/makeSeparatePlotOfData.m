function figHandles = makeSeparatePlotOfData(app,smoothWindow,fignums,usePanels, ...
    varianceType,IQRrange,suppressFigures)
arguments
    app
    smoothWindow = 1;
    fignums = [];
    usePanels=true;
    varianceType = 'STD';  % Plot type standard deviation ('STD') or interquartile range ('IQR'). {'STD', 'IQR'}
    IQRrange = 0.25 % Interquartile range for plotting.
    suppressFigures = false % Hide figures.  Will need to manually make them visible again for viewing.
end

if suppressFigures
    visible = 'off';
    set(0,'DefaultFigureVisible','off')
    if ~isempty(fignums)&&(~iscell(fignums)||~ishandle(fignums{1}))
        disp('Creating new figures -- to reuse figures, fignum must be a cell of figure handles');   
        fignums = [];
    end
else
    visible = 'on';
end

%% This function creates a histogram from the loaded data.
numTimes = length(app.ParEstFitTimesList.Value);
NdModFit = max(1,length(size(app.DataLoadingAndFittingTabOutputs.fitResults.current))-1);
NdDat = length(app.SpeciesForFitPlot.Value);
if NdModFit~=NdDat
    error('Model and Data size do not match in fitting routine')
end

% Plts_to_make = zeros(1,NdModFit);
figHandles = {};
% Density and Cumulative Distributions plots
for DistType = 0:1
    % Choose which figures to plot in
    if isempty(fignums); figHandle = figure('Visible',visible);
    elseif iscell(fignums); figHandle = fignums{DistType+1}; set(0, 'CurrentFigure', figHandle);
    else; figHandle = figure(fignums(DistType+1)); end
    if DistType==0; set(figHandle,'Name','Cumulative Distributions versus Time')
    elseif DistType==1; set(figHandle,'Name','Probability Distributions versus Time')
    end

    % initialize histogram plots
    if numTimes<=4
        subPlotSize = [1,numTimes];
    else
        subPlotSize(2) = ceil(sqrt(numTimes));
        subPlotSize(1) = ceil(numTimes/subPlotSize(2));
    end

    modelHistTime = cell(numTimes,NdDat);
    dataHistTime = cell(numTimes,NdDat);

    % loop over number of time points
    for iTime = 1:numTimes
        % Find index of time in model solution times set
        % iTime = find(strcmp(app.ParEstFitTimesList.Value,app.ParEstFitTimesList.Value{iTime}));
        if usePanels
            subplot(subPlotSize(1),subPlotSize(2),iTime)
        end
        FigHandleHists = gca;

        LegNames = {};
        xm = 0; ym = 0;

        % Switches plotted data, (e.g., x1, x2, x3,...), and legend depending on if the
        % species was selected for plotting.
        for icb=1:NdDat
            cb = max(contains(app.SpeciesForFitPlot.Items,app.SpeciesForFitPlot.Value{icb}));
            % Plts_to_make(cb) = 1;

            indsToSumOver = setdiff([1:NdModFit],icb);
            % soDat = setdiff([1:NdDat],icb);
            speciesName = app.NameTable.Data{icb};

            %% Histograms for the data.
            if cb  % If this species selected for plotting
                dataTensor = double(app.DataLoadingAndFittingTabOutputs.dataTensor);
                if length(size(app.DataLoadingAndFittingTabOutputs.dataTensor))==1
                    H1 = dataTensor;
                else
                    try H1 = squeeze(dataTensor(iTime,:,:,:,:,:,:,:,:,:));
                    catch; whos;
                    end
                end
                if ~isempty(indsToSumOver)
                    H1 = sum(H1,indsToSumOver);
                end
                H1 = H1(:)/sum(H1(:));  % Normalize data before plotting.
                dataHistTime{iTime,icb} = H1;

                % make data histogram plots
                if DistType
                    stairs(FigHandleHists,[0:length(H1)-1],smoothBins(H1,smoothWindow),'linewidth',3);
                    hold(FigHandleHists,'on')
                    ym = max(ym,max(H1));
                else
                    stairs(FigHandleHists,[0:length(H1)-1],cumsum(H1),'linewidth',3);
                    hold(FigHandleHists,'on')
                    ym=1;
                end

                % add name to legend
                LegNames{end+1} = [speciesName,'-data'];
                xm = max(xm,length(H1));

                %% Histograms for the model
                modelTensor = app.DataLoadingAndFittingTabOutputs.fitResults.current;
                inds = modelTensor.subs;
                inds = inds(inds(:,1) == iTime,:);
                vals = modelTensor(inds);
                modelTensor = double(sptensor(inds(:,2:end),vals));
                if ~isempty(indsToSumOver)
                    H1 = squeeze(sum(modelTensor,indsToSumOver));
                else
                    H1 = squeeze(modelTensor);
                end
                modelHistTime{iTime,icb} = H1;

                % make model histogram plots
                if DistType
                    stairs(FigHandleHists,[0:length(H1)-1],smoothBins(H1,smoothWindow),'linewidth',3);
                    hold(FigHandleHists,'on')
                    ym = max(ym,max(H1));
                else
                    stairs(FigHandleHists,[0:length(H1)-1],cumsum(H1),'linewidth',3);
                    hold(FigHandleHists,'on')
                    ym=1;
                end

                % add name to legend
                LegNames{end+1} = [speciesName,'-mod'];
            end

            %% Adjust histograms plot style
            title(['t = ',app.ParEstFitTimesList.Value{iTime}])
            if mod(iTime,subPlotSize(2))==1
                ylabel('Probability')
                %             ylim = get(gca,'ylim');
                %         else
                %             set(gca,'yticklabels',{})
                %             set(gca,'ylim',ylim);
            end
            if iTime>numTimes-subPlotSize(2)
                xlabel('Num. Mol.')
                %         else
                %             set(gca,'xticklabels',{})
            end
            set(gca,'fontsize',15)
        end
    end
    legend(LegNames)
    figHandles{end+1} = figHandle;
end

%% Make trajectory plots for model.
% NdModAll = size(app.NameTable.Data,1);
if isempty(fignums); figHandle = figure('Visible',visible);
elseif iscell(fignums); figHandle = fignums{3}; set(0, 'CurrentFigure', figHandle);
else; figHandle = figure(fignums(3));
end

switch varianceType
    case 'STD'
        set(figHandle,'Name','Mean +/- STD vs. time')
end

meanVarTrajAxis = gca;

tArrayModel = eval(app.FspPrintTimesField.Value);
% if length(tArrayModel)==numTimes
%     tArrayPlot = tArrayModel;
% else
%     tArrayPlot = app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes;
% end

for iTime = 1:length(tArrayModel)
    for j=1:NdDat
        soInds = find(~strcmp(app.SpeciesForFitPlot.Items,app.SpeciesForFitPlot.Value{j}));
        px = app.FspTabOutputs.solutions{iTime}.p;
        if ~isempty(app.FIMTabOutputs.distortionOperator)
            px = app.FIMTabOutputs.distortionOperator.computeObservationDist(px,soInds);
        end
        if isempty(soInds)
            Z = double(px.data);
        else
            Z = double(px.sumOver(soInds).data);
        end
        mnsMod(iTime,j) = [0:length(Z)-1]*Z;
        mns2Mod(iTime,j) = [0:length(Z)-1].^2*Z;

        % Compute the 25% and 75% range from model.
        sumZ = cumsum(Z);
        [~,I] = unique(sumZ);

        if strcmp(varianceType,'IQR')
            if sumZ(1)>IQRrange
                lowIQRmod(iTime,j) = 0;
            else
                lowIQRmod(iTime,j) = interp1(sumZ(I),I-1,IQRrange);
            end
            if sumZ(1)>1-IQRrange
                highIQRmod(iTime,j) = 0;
            else
                highIQRmod(iTime,j) = interp1(sumZ(I),I-1,1-IQRrange);
            end
        end
    end
end
% for iTimeModel = 1:length(T_array)
%      if ~isempty(app.FspTabOutputs.solutions{iTimeModel})
%          for i=1:NdModAll
%              INDS = setdiff([1:NdDat],i);
%
%              % Add effect of PDO.
%              px = app.FspTabOutputs.solutions{iTimeModel}.p;
%              if ~isempty(app.FIMTabOutputs.distortionOperator)
%                  px = app.FIMTabOutputs.distortionOperator.computeObservationDist(px);
%              end
%              if ~isempty(INDS)
%                  mdist{i} = double(px.sumOver(INDS).data);
%              else
%                  mdist{i} = double(px.data);
%              end
%          end
%          for j=1:NdModAll
%              mns(iTimeModel,j) = [0:length(mdist{j})-1]*mdist{j};
%              mns2(iTimeModel,j) = [0:length(mdist{j})-1].^2*mdist{j};
%          end
%      end
%  end
% end
vars = mns2Mod-mnsMod.^2;
cols = ['b','r','g','m','c','k'];
cols2 = [.90 .90  1.00; 1.00 .90 .90; .90 1.00 .90; .60 .60  1.00; 1.00 .60 .60; .60 1.00 .60];
LG = {};

%% TODO - 1:NdDat can be changed to plot species separately (e.g., separate cytGR and nucGR plots)
for iplt=1:NdDat
    % if Plts_to_make(iplt)
    switch varianceType
        case 'STD'
            BD = [mnsMod(:,iplt)'+sqrt(vars(:,iplt)'),mnsMod(end:-1:1,iplt)'-sqrt(vars(end:-1:1,iplt)')];
        case 'IQR'
            BD = [highIQRmod(:,iplt)',lowIQRmod(end:-1:1,iplt)'];

    end
    TT = [tArrayModel(1:end),tArrayModel(end:-1:1)];
    fill(meanVarTrajAxis,TT,BD,cols2(iplt,:));
    hold(meanVarTrajAxis,'on');
    plot(meanVarTrajAxis,tArrayModel,mnsMod(:,iplt),cols(iplt),'linewidth',2);
    LG{end+1} = [char(app.NameTable.Data(iplt)),' Model Mean \pm std'];
    LG{end+1} = [char(app.NameTable.Data(iplt)),' Model Mean'];
    % endc
end
title('Trajectory of means and standard deviations')

%% Add data to trajectory plot.
for iTime = 1:numTimes
    %% TODO - 1:NdDat can be changed to plot species separately (e.g., separate cytGR and nucGR plots)
    for j=1:NdDat
        mnsDat(iTime,j) = [0:length(dataHistTime{iTime,j})-1]*dataHistTime{iTime,j};
        mns2Dat(iTime,j) = [0:length(dataHistTime{iTime,j})-1].^2*dataHistTime{iTime,j};
        % Compute the 25% and 75% range from data.

        sumZ = cumsum(dataHistTime{iTime,j});
        [~,I] = unique(sumZ);

        if sumZ(1)>0.25
            lowIQRmod(iTime,j) = 0;
        else
            lowIQRdat(iTime,j) = interp1(sumZ(I),I-1,0.25);
            highIQRdat(iTime,j) = interp1(sumZ(I),I-1,0.75);
        end

    end
end
varDat = mns2Dat-mnsDat.^2;
T_array = app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes;
%% TODO - 1:NdDat can be changed to plot species separately (e.g., separate cytGR and nucGR plots)
for j=1:NdDat
    switch varianceType
        case 'STD'
            errorbar(meanVarTrajAxis,T_array,...
                mnsDat(:,j),sqrt(varDat(:,j)),[cols(j),...
                'o'],'MarkerSize',12,'MarkerFaceColor',cols(j),'linewidth',3)
        case 'IQR'
            errorbar(meanVarTrajAxis,T_array,...
                mnsDat(:,j),(mnsDat(:,j)-lowIQRdat(:,j)),(highIQRdat(:,j)-mnsDat(:,j)),[cols(j),...
                'o'],'MarkerSize',12,'MarkerFaceColor',cols(j),'linewidth',3)
    end
    LG{end+1} = [char(app.NameTable.Data(j)),' Data mean \pm std'];
end

% Ned = length(size(app.DataLoadingAndFittingTabOutputs.dataTensor))-1;
% for iplt = 1:Ned
% %     if Plts_to_make(iplt)
%         so = setdiff([1:Ned],iplt)+1;
%         dataTensor = double(app.DataLoadingAndFittingTabOutputs.dataTensor);
%         if ~isempty(so)
%             dataTensor = squeeze(sum(dataTensor,so));
%         end
%         P = dataTensor./sum(dataTensor,2);
%         mn = P*[0:size(dataTensor,2)-1]';
%         mn2 = P*([0:size(dataTensor,2)-1]').^2;
%         var = mn2-mn.^2;
%         errorbar(meanVarTrajAxis,app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes,...
%             mn,sqrt(var),[cols(iplt),...
%             'o'],'MarkerSize',12,'MarkerFaceColor',cols(iplt),'linewidth',3)
%         LG{end+1} = [char(app.NameTable.Data(iplt,2)),' Data mean \pm std'];
% %     end
% end
legend(meanVarTrajAxis,LG)

xlabel('Time')
ylabel('Response')
set(meanVarTrajAxis,'fontsize',15)
if max(T_array)>0
    set(meanVarTrajAxis,'xlim',[0,max(T_array)])
end
figHandles{end+1} = figHandle;

%% Make plots of Likelihood functions versus time.
if isempty(fignums);figHandle=figure('Visible',visible);
elseif iscell(fignums); figHandle = fignums{4}; set(0, 'CurrentFigure', figHandle);
else; figHandle=figure(fignums(4));
end
subplot(3,1,1)
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes,app.DataLoadingAndFittingTabOutputs.V_LogLk,'linewidth',2)
hold on
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes,app.DataLoadingAndFittingTabOutputs.perfectMod,'linewidth',2)
% plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,app.DataLoadingAndFittingTabOutputs.perfectModSmoothed,'linewidth',2)
ylabel('Log(L(D(t)|M)')
% legend({'log(L) - best fit for model','log(L) - theoretical limit','log(L) - theoretical limit (smoothed)'})
legend({'log(L) - best fit for model','log(L) - theoretical limit'})

subplot(3,1,2)
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes,app.DataLoadingAndFittingTabOutputs.numCells,'linewidth',2)
ylabel('# Cells')

subplot(3,1,3)
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes,...
    (app.DataLoadingAndFittingTabOutputs.V_LogLk - app.DataLoadingAndFittingTabOutputs.perfectMod)./...
    app.DataLoadingAndFittingTabOutputs.numCells,'linewidth',2)

ylabel('\Delta Log(L) / Cell')
xlabel('Time')
for i=1:3
    subplot(3,1,i)
    set(gca,'fontsize',15);
end
figHandles{end+1} = figHandle;

%% Make Joint Density Plots if 2 or more species in data set
% if NdDat==2
%     ModelTensor = app.DataLoadingAndFittingTabOutputs.fitResults.current;
%     DataTensor = app.DataLoadingAndFittingTabOutputs.fitResults.currentData;
%     figHandle = figure;
%     for iTime = 1:size(ModelTensor,1)
%         Z = log10(squeeze(ModelTensor(iTime,:,:))+1e-6);
%         subplot(subPlotSize(1),subPlotSize(2),iTime)
%         contourf(Z)
%     end
%     figHandle = figure;
%     for iTime = 1:size(ModelTensor,1)
%         Z = log10(squeeze(DataTensor(iTime,:,:))+1e-6);
%         subplot(subPlotSize(1),subPlotSize(2),iTime)
%         contourf(Z)
%     end
% end
end

function sb = smoothBins(x,bnsz)
sb = 0*x;
for i=1:length(x)
    j = bnsz*floor(i/bnsz);
    sb(j+1:min(j+bnsz,length(x))) = sb(j+1:min(j+bnsz,length(x)))+x(i);
end
end
