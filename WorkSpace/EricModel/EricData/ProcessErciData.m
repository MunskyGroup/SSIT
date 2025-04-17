Orig = readtable('Complete_dataframe_Ron_010224.csv');
X = Orig;
%%
X.Properties.VariableNames

%% Estimate Cytoplasm Size for Un-segmented cytoplasms.
figure(1); clf;
J = isnan(X.Cyto_Area);
CytoAreas = X.Cyto_Area(~isnan(X.Cyto_Area));
NucAreas = X.Nuc_Area(~isnan(X.Cyto_Area));
scatter(NucAreas,CytoAreas,'bx')
P = polyfit(NucAreas,CytoAreas,2);
x = linspace(0,max(NucAreas),100)';
y = polyval(P,x);
X.Cyto_Area(isnan(X.Cyto_Area)) = polyval(P,X.Nuc_Area(isnan(X.Cyto_Area)));

hold on
scatter(X.Nuc_Area(J),X.Cyto_Area(J),'xr')
plot(x,y,'k--','LineWidth',3)

%% Gate on 25% to 75% of nucleus size.
NA = X.Nuc_Area;
sortNA = sort(NA);
nucLow = sortNA(floor(0.25*length(sortNA)));
nucHigh = sortNA(ceil(0.75*length(sortNA)));
J = X.Nuc_Area>=nucLow&X.Nuc_Area<=nucHigh;
X = X(J,:);

%% Estimate volume ratio between nucleus and cytoplasm (using segemented cells)
sortNA = sort(NucAreas);
nucLow = sortNA(floor(0.25*length(sortNA)));
nucHigh = sortNA(ceil(0.75*length(sortNA)));
J = NucAreas>=nucLow&NucAreas<=nucHigh;
ratio = mean(NucAreas(J))/mean(CytoAreas(J))

%% Compute "normalized" GR in nuc and cytoplasm.
figure(2); clf;
J = contains(X.Condition,'GR_timesweep');
binsCyt = 15; 
binsNuc = 30;

threshold = 0.005;
GRDat = X(J,:);


nGRcells = size(GRDat,1);

sortGRnuc = sort(GRDat.Nuc_GR_avg_int);
nucGR95 = sortGRnuc(ceil((1-threshold)*length(sortGRnuc)));
nucGR05 = sortGRnuc(ceil(threshold*length(sortGRnuc)));
nucGRnorm = round(binsNuc*max(0,min(1,(GRDat.Nuc_GR_avg_int-nucGR05)/(nucGR95-nucGR05))));
subplot(1,3,1);histogram(round(nucGRnorm),'BinEdges',[0:binsNuc])
GRDat.normgrnuc = nucGRnorm;
title('Nuclear Distribution')

sortGRcyt = sort(GRDat.Cyto_GR_avg_int);
cytGR95 = sortGRcyt(ceil((1-threshold)*length(sortGRcyt)));
cytGR05 = sortGRcyt(ceil(threshold*length(sortGRcyt)));
cytGRnorm = round(binsCyt*max(0,min(1,(GRDat.Cyto_GR_avg_int-cytGR05)/(cytGR95-cytGR05))));
subplot(1,3,2);histogram(round(cytGRnorm),'BinEdges',[0:binsCyt])
GRDat.normgrcyt = cytGRnorm;
title('Cytoplasm Distribution')

%% Randomly assign zero minute timepoints to other GR data.
jZero = find(GRDat.Dex_Conc == 0);
iZero = ceil(3*rand(size(jZero)));
dexConc = [1,10,100];
GRDat.Dex_Conc(jZero) = dexConc(iZero);

%%
    
jointHist = zeros(binsNuc+1,binsCyt+1);
for i=1:length(nucGRnorm)
    jointHist(nucGRnorm(i)+1,cytGRnorm(i)+1) = jointHist(nucGRnorm(i)+1,cytGRnorm(i)+1) + 1;
end
subplot(1,3,3);contourf(jointHist)
ylabel('Nuclear')
xlabel('Cytoplasm')

GRDatReduced = array2table([GRDat.Cell_id,GRDat.Time_index,GRDat.Dex_Conc,GRDat.normgrnuc,GRDat.normgrcyt],...
'variableNames',{'Cell_id','time','Dex_Conc','normgrnuc','normgrcyt'});

writetable(GRDatReduced,'Gated_dataframe_Ron_020224_NormalizedGR_bins.csv')

return

% %% Extract Data for GR experiments
% repStrings = unique(GRDat.Replica);
% times = unique(GRDat.Time_index);
% concs = [0,1,10,100];
% 
% for iconc = 1:4
%     for irep = 1:length(repStrings)
%         for itime = 1:length(times)
%             GRDatCell{iconc,irep,itime} = GRDat(([GRDat.Replica{:}]'==repStrings{irep})&...
%                 (GRDat.Time_index==times(itime))...
%                 &(GRDat.Dex_Conc==concs(iconc)),:);
% 
%             medianNuc(iconc,irep,itime) = median(GRDatCell{iconc,irep,itime}.Nuc_GR_avg_int);
%             medianCyt(iconc,irep,itime) = median(GRDatCell{iconc,irep,itime}.Cyto_GR_avg_int);
% 
%         end
%     end
% end
% 
% %% Compute Total GR in Nucleus and Cytoplasm and bin.
% for iconc = 1:4
%     for irep = 1:length(repStrings)
%         for itime = 1:length(times)
%             GRDatCell{iconc,irep,itime}
% 
%         end
%     end
% end
% %%
% edgesNuc = linspace(0,8000,50);
% edgesCyt = linspace(0,4000,50);
% for iconc = 1:4
% 
%     figure(10*iconc+1);clf;
%     figure(10*iconc+2);clf;
%     for irep = 1:length(repStrings)
%         for itime = 1:length(times)
% 
%             if itime==1
%                 DT = GRDatCell{1,irep,itime};
%             else
%                 DT = GRDatCell{iconc,irep,itime};
%             end
% 
%             figure(10*iconc+1);
%             subplot(3,4,itime)
%             histogram(DT.Nuc_GR_avg_int,edgesNuc,'Normalization','probability')
%             set(gca,'xlim',[0,8000],'ylim',[0,0.2])
%             title([num2str(times(itime)),' min'])
%             hold on
%             plot(mean(medianNuc(1,:,1))*[1,1],[0,0.2],'k--','LineWidth',2)
% 
%             figure(10*iconc+2);
%             subplot(3,4,itime)
%             histogram(DT.Cyto_GR_avg_int,edgesCyt,'Normalization','probability')
%             set(gca,'xlim',[0,4000],'ylim',[0,0.3])
%             hold on
%             plot(mean(medianCyt(1,:,1))*[1,1],[0,0.3],'k--','LineWidth',2)
%             title([num2str(times(itime)),' min'])
%         end
% 
%     end
% end
% 
% 
% % for t = 
% 
% 
