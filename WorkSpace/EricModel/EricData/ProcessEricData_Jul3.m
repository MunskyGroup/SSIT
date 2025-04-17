%% Load original data
Orig = readtable('Complete_dataframe_Ron_010224.csv');
X = Orig;
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

%% Save Copy of all data after nuc size gating.
writetable(X,'Data010224_Gated_On_Nuc.csv')

%% Estimate volume ratio between nucleus and cytoplasm (using segemented cells)
sortNA = sort(NucAreas);
nucLow = sortNA(floor(0.25*length(sortNA)));
nucHigh = sortNA(ceil(0.75*length(sortNA)));
J = NucAreas>=nucLow&NucAreas<=nucHigh;
ratio = mean(NucAreas(J))/mean(CytoAreas(J))

%% Compute "normalized" GR concentrations in nuc and cytoplasm.
figure(2); clf;
J = contains(X.Condition,'GR_timesweep');
GRDat = X(J,:);

%% Remove timepoints with only 1 replica
J = (GRDat.Time_index~=20)&(GRDat.Time_index~=40)&(GRDat.Time_index~=60)&(GRDat.Time_index~=90)&(GRDat.Time_index~=120);
GRDat = GRDat(J,:)

% Specify the number of bins for cytoplasm concentration.
binsCyt = 15; 

% Threshold for saturation
threshold = 0.01;

% % Use t=0 data to define binning for cytoplasmic GR data
% GRDat_t0 = GRDat(GRDat.Time_index==0,:);

% compute total GR.
totalGR = GRDat.Cyto_GR_avg_int + ratio*GRDat.Nuc_GR_avg_int;

% Bin the cytopasmic GR data
sortGRcyt = sort(GRDat.Cyto_GR_avg_int);
maxCytGR = sortGRcyt(ceil((1-threshold)*length(sortGRcyt)));
minThreshold = sortGRcyt(ceil(threshold*length(sortGRcyt)));
cytGRnorm = round(binsCyt*max(0,min(1,(GRDat.Cyto_GR_avg_int-minThreshold)/(maxCytGR-minThreshold))));
subplot(1,3,2);histogram(round(cytGRnorm),'BinEdges',[0:binsCyt])
GRDat.normgrcyt = cytGRnorm;
title('Cytoplasm Distribution')

% Compute corresponding levels of nuclear GR
% First scale by volume to get the total amount in nucleus
GRnucLevel = GRDat.Nuc_GR_avg_int*ratio;
nucGRnorm = round(binsCyt*max(0,min(1,(GRnucLevel-minThreshold)/(maxCytGR-minThreshold))));
subplot(1,3,1);histogram(round(nucGRnorm),'BinEdges',[0:max(round(nucGRnorm))])
GRDat.normgrnuc = nucGRnorm;
title('Nuclear Distribution')

%% Randomly assign zero minute timepoints to other GR data.
jZero = find(GRDat.Dex_Conc == 0);
iZero = ceil(3*rand(size(jZero)));
dexConc = [1,10,100];
GRDat.Dex_Conc(jZero) = dexConc(iZero);

%% Plot and save results.
jointHist = zeros(max(round(nucGRnorm))+1,max(round(cytGRnorm))+1);
for i=1:length(nucGRnorm)
    jointHist(nucGRnorm(i)+1,cytGRnorm(i)+1) = jointHist(nucGRnorm(i)+1,cytGRnorm(i)+1) + 1;
end
subplot(1,3,3);contourf(jointHist)
ylabel('Nuclear')
xlabel('Cytoplasm')

GRDatReduced = array2table([GRDat.Cell_id,GRDat.Time_index,GRDat.Dex_Conc,GRDat.normgrnuc,GRDat.normgrcyt],...
'variableNames',{'Cell_id','time','Dex_Conc','normgrnuc','normgrcyt'});

writetable(GRDatReduced,'Gated_dataframe_Ron_020224_NormalizedGR_bins.csv')

