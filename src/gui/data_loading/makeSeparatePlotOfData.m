function makeSeparatePlotOfData(app,smoothWindow,fignums,usePanels)
arguments
    app
    smoothWindow = 5;
    fignums = [];
    usePanels=true;
end

%% This function creates a histogram from the loaded data.
NT = length(app.ParEstFitTimesList.Value);
% NdMod = max(1,length(size(app.DataLoadingAndFittingTabOutputs.fitResults.current))-1);
NdMod = size(app.NameTable.Data,1);
NdDat = length(app.SpeciesForFitPlot.Value);

Plts_to_make = zeros(1,NdMod);

for DistType = 0:1
    if isempty(fignums)
        figure
    else
        figure(fignums(DistType+1))
    end
    if NT<=4
        subPlotSize = [1,NT];
    else
        subPlotSize(2) = ceil(sqrt(NT));
        subPlotSize(1) = ceil(NT/subPlotSize(2));
    end

    
    for iTime = 1:NT
        it = find(strcmp(app.ParEstFitTimesList.Value,app.ParEstFitTimesList.Value{iTime}));
        if usePanels
            subplot(subPlotSize(1),subPlotSize(2),iTime)
        end
        FNHists = gca;
%         hold off
        
        L = {};
        xm = 0; ym = 0;
        
        % Switches plotted data, (x1, x2, x3,...), and legend depending on if the
        % species was selected for plotting.
        for icb=1:NdDat
            cb = contains(app.SpeciesForFitPlot.Items,app.SpeciesForFitPlot.Value{icb});
            Plts_to_make(cb) = 1;
            
            soMod = setdiff([1:NdMod],find(cb~=0));
            soDat = setdiff([1:NdDat],icb);
            le = app.NameTable.Data{icb,2};
            
%             if cb
                matTensor = double(app.DataLoadingAndFittingTabOutputs.dataTensor);
                if length(size(app.DataLoadingAndFittingTabOutputs.dataTensor))==1
                    H1 = matTensor;
                else
                    try
                        H1 = squeeze(matTensor(it,:,:,:,:,:,:,:,:,:));
                    catch
                        whos
                    end
                end
                if ~isempty(soDat)
                    H1 = sum(H1,soDat);
                end
                H1 = H1(:)/sum(H1(:));  % Normalize data before plotting.
                
                if DistType
%                     stairs(FNHists,[0:length(H1)-1],smoothBins(H1,smoothWindow)); hold on
                    stairs(FNHists,[0:length(H1)-1],smoothBins(H1,smoothWindow),'linewidth',3);
                    hold(FNHists,'on')
                    ym = max(ym,max(H1));
                else
                    stairs(FNHists,[0:length(H1)-1],cumsum(H1),'linewidth',3);
                    hold(FNHists,'on')
                    ym=1;
                end
                L{end+1} = [le,'-data'];
                xm = max(xm,length(H1));
                
                if ~isempty(soMod)
                    M = squeeze(app.DataLoadingAndFittingTabOutputs.fitResults.current(it,:,:,:,:,:,:,:,:,:,:,:,:));
                    H1 = squeeze(sum(M,soMod));
                else
                    H1 = squeeze(app.DataLoadingAndFittingTabOutputs.fitResults.current(it,:,:,:,:,:,:,:,:,:,:,:,:));
                end
                % H1(end+1)=0;
                if DistType
                    stairs(FNHists,[0:length(H1)-1],smoothBins(H1,smoothWindow),'linewidth',3);
                    hold(FNHists,'on')
                    ym = max(ym,max(H1));
                else
                    stairs(FNHists,[0:length(H1)-1],cumsum(H1),'linewidth',3);
                    hold(FNHists,'on')
                    ym=1;
                end
                L{end+1} = [le,'-mod'];
%             end
        end
        title(['t = ',app.ParEstFitTimesList.Value{iTime}])
%         if mod(iTime,subPlotSize(2))==1
            ylabel('Probability')
%             ylim = get(gca,'ylim');
%         else
%             set(gca,'yticklabels',{})
%             set(gca,'ylim',ylim);
%         end
%         if iTime>NT-subPlotSize(2)
            xlabel('Mol. Count')
%         else
%             set(gca,'xticklabels',{})
%         end
set(gca,'fontsize',15)
    end
end


%% Make a trajectory plot for model.
NdMod = size(app.NameTable.Data,1);
if isempty(fignums)
    figure
else
    figure(fignums(3))
end
FNTraj = gca;
T_array = eval(app.FspPrintTimesField.Value);
for it = 1:length(T_array)
    if ~isempty(app.FspTabOutputs.solutions{it})
        for i=1:NdMod
            INDS = setdiff([1:NdMod],i);

            % Add effect of PDO.
            px = app.FspTabOutputs.solutions{it}.p;
            if ~isempty(app.FIMTabOutputs.distortionOperator)
                px = app.FIMTabOutputs.distortionOperator.computeObservationDist(px);
            end
            if ~isempty(INDS)
                mdist{i} = double(px.sumOver(INDS).data);
            else
                mdist{i} = double(px.data);
            end
        end
        for j=1:NdMod
            mns(it,j) = [0:length(mdist{j})-1]*mdist{j};
            mns2(it,j) = [0:length(mdist{j})-1].^2*mdist{j};
        end
    end
end
vars = mns2-mns.^2;
cols = ['b','r','g','m','c','k'];
cols2 = [.90 .90  1.00; 1.00 .90 .90; .90 1.00 .90; .60 .60  1.00; 1.00 .60 .60; .60 1.00 .60];
LG = {};
for iplt=1:NdMod
    if Plts_to_make(iplt)
        BD = [mns(:,iplt)'+sqrt(vars(:,iplt)'),mns(end:-1:1,iplt)'-sqrt(vars(end:-1:1,iplt)')];
        TT = [T_array(1:end),T_array(end:-1:1)];
        fill(FNTraj,TT,BD,cols2(iplt,:));
        hold(FNTraj,'on');
        plot(FNTraj,T_array,mns(:,iplt),cols(iplt),'linewidth',2);
        LG{end+1} = [char(app.NameTable.Data(iplt,2)),' Model Mean \pm std'];
        LG{end+1} = [char(app.NameTable.Data(iplt,2)),' Model Mean'];
    end
end
title('Trajectory of means and standard deviations')

%% Add data to trajectory plot.
Ned = length(size(app.DataLoadingAndFittingTabOutputs.dataTensor))-1;
for iplt = 1:Ned
%     if Plts_to_make(iplt)
        so = setdiff([1:Ned],iplt)+1;
        matTensor = double(app.DataLoadingAndFittingTabOutputs.dataTensor);
        if ~isempty(so)
            matTensor = squeeze(sum(matTensor,so));
        end
        P = matTensor./sum(matTensor,2);
        mn = P*[0:size(matTensor,2)-1]';
        mn2 = P*([0:size(matTensor,2)-1]').^2;
        var = mn2-mn.^2;
%         plot(FNTraj,app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes,mn,[cols(iplt),...
%             'o'],'MarkerSize',12,'MarkerFaceColor',cols(iplt))
        errorbar(FNTraj,app.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes,...
            mn,sqrt(var),[cols(iplt),...
            'o'],'MarkerSize',12,'MarkerFaceColor',cols(iplt),'linewidth',3)
        LG{end+1} = [char(app.NameTable.Data(iplt,2)),' Data mean \pm std'];
%     end
end
legend(FNTraj,LG)

xlabel('Time')
ylabel('Response')
set(FNTraj,'fontsize',15)
if max(T_array)>0
    set(FNTraj,'xlim',[0,max(T_array)])
end
%%

if isempty(fignums)
    figure
else
    figure(fignums(4))
end
subplot(3,1,1)
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,app.DataLoadingAndFittingTabOutputs.V_LogLk,'linewidth',2)
hold on
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,app.DataLoadingAndFittingTabOutputs.perfectMod,'linewidth',2)
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,app.DataLoadingAndFittingTabOutputs.perfectModSmoothed,'linewidth',2)
ylabel('Log(L(D(t)|M)')
legend({'log(L) - best fit for model','log(L) - theoretical limit','log(L) - theoretical limit (smoothed)'})

subplot(3,1,2)
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,app.DataLoadingAndFittingTabOutputs.numCells,'linewidth',2)
ylabel('# Cells')

subplot(3,1,3)
plot(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times,...
    (app.DataLoadingAndFittingTabOutputs.V_LogLk - app.DataLoadingAndFittingTabOutputs.perfectMod)./...
    app.DataLoadingAndFittingTabOutputs.numCells,'linewidth',2)

ylabel('\Delta Log(L) / Cell')
xlabel('Time')
for i=1:3
    subplot(3,1,i)
    set(gca,'fontsize',15);
end


end

function sb = smoothBins(x,bnsz)
sb = 0*x;
for i=1:length(x)
    j = bnsz*floor(i/bnsz);
    sb(j+1:min(j+bnsz,length(x))) = sb(j+1:min(j+bnsz,length(x)))+x(i);
end
end
