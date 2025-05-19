%% Simulate diffusion
clear
close all 
addpath(genpath('../../../src'));
addpath("RNA_Modeling")

%% Generate Simulated data
run_sim = true;

kon = 1e-3;
koff = 7.5e-5;
w = 0.0025;
kex = 750;
kr = 7.5e4;
D = [0.01,5,4];
gam =[0.035;0.0025;0.001];
Ncells = 100;
a = 53*ones(Ncells,1);
b = 39*ones(Ncells,1);
u = 2*pi*rand(Ncells,1);
posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data
makePlot = 1;
fileName = 'demo1.csv';

if run_sim
    stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)
end

[intensTS_Sims,binCounts_Sims,~,~,~,binTS_Sims,binCounts_SimsDim,binCounts_SimsBright,...
    hisCountsBright_Sims,hisCountsDim_Sims,XcRot_Sims,iTS_Sims,XCVcRot_Sims,XYfit_Sims] = ...
    processSpotPositionData(fileName,[],a,b,0);


%% Simple Model Fit

figure(1)
plot(XcRot_Sims(:,1),XcRot_Sims(:,2),'o'); hold on
plot(XCVcRot_Sims(:,1),XCVcRot_Sims(:,2),'k--','lineWidth',1.5); hold on
% plot(XcRot_Sims(iTS_Sims,1),XcRot_Sims(iTS_Sims,2),'rx','MarkerSize',12,'lineWidth',2)
plot(-posnTS(18,1),-posnTS(18,2),'rx','MarkerSize',12,'lineWidth',2)


figure(2);
edges = [0:20:600];
spotsPerCellMod = sum(binCounts_Sims);
spotsPerCellMod(isnan(spotsPerCellMod))=0;
Hmodel = histcounts(spotsPerCellMod,edges,'Normalization','pdf');
plot(edges(1:end-1)+10,Hmodel);
legend('Data','Model')



for icase = 2
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


%% Save Bin Counts as CSV
binned_table = array2table(binCounts_Sims', ...
    'VariableNames', {'bin_1', 'bin_2', 'bin_3', 'bin_4', 'bin_5', 'bin_6', 'bin_7', 'bin_8', 'bin_9', 'bin_10', 'bin_11'});

writetable(binned_table, 'sim.csv')








