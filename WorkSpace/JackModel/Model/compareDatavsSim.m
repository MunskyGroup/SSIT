clear all

%% Process Experimental data
dataFile = 'dataframe_MS2-CY5_Cyto543_560_woStim.csv';
maskFile = 'polygons_wo.mat';
[intensTS_Data,binCounts_Data,abEllipse_Data,areaEllipse_Data,...
    ecentricity_Data,binTS_Data,binCounts_DataDim,binCounts_DataBright,...
    hisCountsBright_Data,hisCountsDim_Data,XcRot_Data,iTS_Data,XCVcRot_Data,XYfit_Data] = ...
    processSpotPositionData(dataFile,maskFile,[],[],0);

%% Generate Simulated data
kon = 1e-3;
koff = 7.5e-5;
w = 0.0025;
kex = 750;
kr = 7.5e4;
D = [0.01,5,4];
gam =[0.035;0.0025;0.001];
Ncells = 2000;
a = 53*ones(Ncells,1);
b = 39*ones(Ncells,1);
u = 2*pi*rand(Ncells,1);
posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data
% makePlot = 1;`!!qSim.csv';

stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

%% Process Simulated Data.
[intensTS_Sims,binCounts_Sims,~,~,~,binTS_Sims,binCounts_SimsDim,binCounts_SimsBright,...
    hisCountsBright_Sims,hisCountsDim_Sims,XcRot_Sims,iTS_Sims,XCVcRot_Sims,XYfit_Sims] = ...
    processSpotPositionData(fileName,[],a,b,0);

%% Scatter Plots
figure(1)
subplot(2,1,1)
plot(XcRot_Data(:,1),XcRot_Data(:,2),'o'); hold on
plot(XCVcRot_Data(:,1),XCVcRot_Data(:,2),'k--','lineWidth',1.5); hold on
plot(XcRot_Data(iTS_Data,1),XcRot_Data(iTS_Data,2),'rx','MarkerSize',12,'lineWidth',2)
title('Experimental Data')

subplot(2,1,2)
plot(XcRot_Sims(:,1),XcRot_Sims(:,2),'o'); hold on
plot(XCVcRot_Sims(:,1),XCVcRot_Sims(:,2),'k--','lineWidth',1.5); hold on
% plot(XcRot_Sims(iTS_Sims,1),XcRot_Sims(iTS_Sims,2),'rx','MarkerSize',12,'lineWidth',2)
plot(-posnTS(18,1),-posnTS(18,2),'rx','MarkerSize',12,'lineWidth',2)

title('Simulated Data')
% iCell+1,:
%  for i=1:15
%     binCounts(i,iCell+1) = sum((XcRot(:,1)>=binEdges(i)&XcRot(:,1)<=binEdges(i+1)));
%  end

%% Make figures to compare model and data.
% Histograms of mRNA counts
figure(2); clf
edges = [0:20:600];
spotsPerCellData = sum(binCounts_Data);
spotsPerCellData(isnan(spotsPerCellData))=0;
Hdata = histcounts(spotsPerCellData,edges,'Normalization','pdf');
spotsPerCellMod = sum(binCounts_Sims);
spotsPerCellMod(isnan(spotsPerCellMod))=0;
Hmodel = histcounts(spotsPerCellMod,edges,'Normalization','pdf');
plot(edges(1:end-1)+10,Hdata,edges(1:end-1)+10,Hmodel);
legend('Data','Model')

%% Make Figure to compare distribution for Number of Nascent mRNA per TS. % what is this showing
figure
maxNascent = max(intensTS_Data);
nascentHistData = histcounts(intensTS_Data,[0:maxNascent+1]);
nascentHistData=nascentHistData/sum(nascentHistData);
nascentHistSims = histcounts(intensTS_Sims,[0:maxNascent+1]);
nascentHistSims=nascentHistSims/sum(nascentHistSims);
plot([0:maxNascent],cumsum(nascentHistData),[0:maxNascent],cumsum(nascentHistSims),'--','LineWidth',3);
legend('Data','Model')
title('Compare distribution of nascent mRNA')
KS_Nascent_Dist = max(abs(cumsum(nascentHistData)-cumsum(nascentHistSims)))

%% Make Figure to compare distribution for Number of Mature mRNA per TS.
figure
maxMature = max(sum(binCounts_Data));
matureHistData = histcounts(sum(binCounts_Data),[0:maxMature+1]);
matureHistData=matureHistData/sum(matureHistData);
matureHistSims = histcounts(sum(binCounts_Sims),[0:maxMature+1]);
matureHistSims=matureHistSims/sum(matureHistSims);
plot([0:maxMature],cumsum(matureHistData),[0:maxMature],cumsum(matureHistSims),'--','LineWidth',3);
legend('Data','Model')
title('Compare distribution of mature mRNA')
KS_Mature_Dist = max(abs(cumsum(matureHistData)-cumsum(matureHistSims)))

%% Histograms of Spots in Bins along major axis. % How do these change as a function of sim params
% Spots in spatial bins.
figure; clf

for icase = 1:2 
    switch icase
        case 1
            binCounts=binCounts_Data;
            binTS=binTS_Data;
            lsp = '-';
        case 2
            binCounts=binCounts_Sims;
            binTS=binTS_Sims;
            lsp = '--';
        case 3
            binCounts=binCounts_SimsDim;
            binTS=binTS_Sims;
            lsp = '--';
        case 4
            binCounts=binCounts_SimsBright;
            binTS=binTS_Sims;
            lsp = '--';
        case 5
            binCounts=binCounts_DataDim;
            binTS=binTS_Data;
            lsp = '-';
        case 6
            binCounts=binCounts_DataBright;
            binTS=binTS_Data;
            lsp = '-';
    end

    for ibinloc = 1:5
        J = (binTS == ibinloc);
        binCountsLoc{ibinloc} = binCounts(:,J);
        N = size(binCountsLoc{ibinloc},2);
        nBins = size(binCountsLoc{ibinloc},1);
        subplot(3,2,ibinloc)
        spatialConcentration{icase} = sum(binCountsLoc{ibinloc},2)/N;
        stairs([0:nBins-1],spatialConcentration{icase},lsp,'lineWidth',3);

        hold on
%       title(['RNA Counts vs. Bin Location Cell 18'],'FontSize',20)
        title('TS Bin Location = ',num2str(ibinloc),'FontSize',20)
%       plot([ibinloc-0.5,ibinloc+0.5],[0,0],'rx','MarkerSize',12,'lineWidth',2)
        plot([ibinloc],[36.6],'rx','MarkerSize',12,'lineWidth',2)
        xlabel('Spatial Position','FontSize',14); ylabel('Number of RNA','FontSize',14)
        
    end
end

%             legend for Fig 6
%             h = legend('RNA','TS');
%             set(h,'FontSize',15);
%     
%             legend for Fig 6
%             h = legend('Total Spots','TS');
%             set(h,'FontSize',15);
%         
%            legend for Fig 7
%             h = legend('Total Spots','TS','Dim Spots','','Bright Spots');
%             set(h,'FontSize',15);
%     
%             legend for Fig 8
                h = legend('Experimental Data','TS','Simulated Data');
                set(h,'FontSize',15);
%             totalError = sum((spatialConcentration{1} - spatialConcentration{2}).^2)
%     
%             legend for Fig 9
%             h = legend('Total Spots Data','TS','Total Spots Sim','','Dim Spots Sim','','Bright Spots Sim','','Dim Spots Data','','Bright Spots Data');
%             set(h,'FontSize',15);

%% Make Figures to compare distance of bright and Dim mRNA from TS % this needs to be completely revisited
figure
hisCountsBrightSims = mean(hisCountsBright_Sims,'omitnan');
hisCountsBrightSims = hisCountsBrightSims/sum(hisCountsBrightSims);
hisCountsDimSims = mean(hisCountsDim_Sims,'omitnan');
hisCountsDimSims = hisCountsDimSims/sum(hisCountsDimSims);
hisCountsBrightDat = mean(hisCountsBright_Data,'omitnan');
hisCountsBrightDat = hisCountsBrightDat/sum(hisCountsBrightDat);
hisCountsDimDat = mean(hisCountsDim_Data,'omitnan');
hisCountsDimDat = hisCountsDimDat/sum(hisCountsDimDat);

subplot(1,2,1)
plot(hisCountsBrightDat,'r','LineWidth',3);hold on;
plot(hisCountsDimDat,'b','LineWidth',3);hold on;
plot(hisCountsBrightSims,'r--','LineWidth',3);hold on;
plot(hisCountsDimSims,'b--','LineWidth',3);hold on;
legend('Data - Bright','Data - Dim','Model - Bright','Model - Dim')
title('Distance Between TS and RNAs')
xlabel('r')

subplot(1,2,2)
plot(cumsum(hisCountsBrightDat),'r','LineWidth',3);hold on;
plot(cumsum(hisCountsDimDat),'b','LineWidth',3);hold on;
plot(cumsum(hisCountsBrightSims),'r--','LineWidth',3);hold on;
plot(cumsum(hisCountsDimSims),'b--','LineWidth',3);hold on;
legend('Data - Bright','Data - Dim','Model - Bright','Model - Dim')
title('Distance Between TS and RNAs')
xlabel('r')


%% Make histogram of spot intensities in the data
figure
Xdat = importdata(dataFile);
j = find(strcmp(Xdat.colheaders,'spot_int_ch_3'));
histogram(Xdat.data(:,j),linspace(-2000,2000,100))