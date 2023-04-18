%% Vo_et_all_PlottingResults
% This script provides all codes needed to generate the figures for the
% fitting and FIM analysis for the HIV promoter example in for Vo et al,
% 2023.  To generate the figures, you will need to run the individual cells
% one at a time.  

%% Preliminary Setup to select model and specifiy setting used in the fitting.
clear all
addpath('../')
TMP = SSIT(); clear TMP % Load definition of the SSIT class for later use.
FN = 'Parameter_Fit_Files/FISHTrue_2StateBurst_4Pars_Prior2x';
vars.doFit = 0;
vars.pdoTimes = [0,300];
vars.priorScale = 2;
vars.muLog10Prior = [-4,-4,log10(0.2),log10(7.12),log10(0.012)]';
vars.sigLog10Prior = [1,1,.5,.5,.5]'*vars.priorScale;
vars.makePlots = false;
vars.modelVarsToFit = [2:5];
vars.initialFspBounds = [0 0 0 0 3 3 3 1200];
vars.timeSet = [0,18,300];
vars.modelChoice = '2stateBurst';
vars.covWithPrior = true;

%% Make Figure to compare PDO statistics and select PDO (Figs 6,7)
NCells = [135 96 62];
kPars = [3,4,7]';
makeFiguresForPDOs(vars,FN,NCells,kPars);

%% Make Initial Fit Figures (Figs 7B. Left 4 Columns)
close all
fitResults = makeFitFigures(vars,[FN,'_0_300']);

%% Make Final Fit Figures (Not shown in paper, but may be useful)
close all
makeFitFigures(vars,[FN,'_0_18_300']);

%% Make Initial MH Scatter Plots (Figs 7B. Right 2 Columns)
close all
makeMHandFimPlots(vars,FN,[FN,'_0_300'])

%% Make Final MH Scatter Plots (Not shown in paper, but may be useful)
close all
makeMHandFimPlots(vars,FN,[FN,'_0_18_300'])

%% Make Tex Table for Initial parameters (Table 3)
makeTexParTable([FN,'_0_300'])

%% Make Tex Table for Final parameters (Table 4)
makeTexParTable([FN,'_0_18_300'])

%% Make Plot of FIM vs time of third experiment (Fig 8a)
vars.covWithPrior = true;
vars.FIMShowErrorBars = true;
vars.FIMnCells = [135,0,62];
makeFimExptDesign(vars,FN,[-10,-7;-8,-5;-7,-4;-6,-1]-3)

%% Compare FIM for different experiment designs (Figs 8b,c,d)
close all
NCells = [135 96 62];
[detCovFIM,detCovMH,FIM_0_300] = makeFimExpComparisons(FN,NCells,vars);
log10(detCovFIM)
log10(detCovMH)

%% Make Fit Figures for Simplified Model (Figure S11, Left 4 Columns)
vars.modelVarsToFit = [2:4];
vars.initialFspBounds = [0 0 0 3 3 1200];
vars.muLog10Prior=[-4,-4,log10(0.2),log10(0.012)]';
vars.sigLog10Prior = [1,1,0.5,.5]'*vars.priorScale;
vars.modelChoice = '2statePoisson';
FN = 'Parameter_Fit_Files/FISHTrue_3Pars_TwoStatePoisson';

close all
fitResultsSimp = makeFitFigures(vars,[FN,'_0_300']);
fitMLE = sum(fitResultsSimp.fitLikelihoods(:,[1,3]),2);
fitMLE([1,3,4,5])
BIC = -2*fitMLE([1,3,4,5])+4*log(135+62)
predMLE = sum(fitResultsSimp.fitLikelihoods(:,[2]),2);
predMLE([1,3,4,5,7,8,9])

%% Make MH plots for simplified model (Figure S11, Right 2 Columns)
close all
ax1 = [-4.2,-3.2; -1,1; -3.75,-1.75];
ax4 = [-5.5,-3; -1,2; -3.5,1];
makeMHandFimPlots(vars,FN,[FN,'_0_300'],true,ax1,ax4)

%% Make Plot of FIM vs time of third experiment (Fig S12a)
vars.covWithPrior = true;
vars.FIMShowErrorBars = true;
vars.FIMnCells = [135,0,62];
makeFimExptDesign(vars,FN,[-11,-9;-10,-7;-10,-7;-8,-4])

%% Compare FIM for different experiment designs (Figs S12b,c,d)
close all
NCells = [135 96 62];
[detCovFIM,detCovMH] = makeFimExpComparisons(FN,NCells,vars);
log10(detCovFIM)
log10(detCovMH)

%% Functions
function [Final_Total,Sigs_Total,Final_STD,FIM_0_300] = evaluateFIMDesigns(FIMR,iPDO,NCells,vars)
arguments
    FIMR
    iPDO
    NCells
    vars
end

covPrior0 = vars.sigLog10Prior.^2;
pars = vars.modelVarsToFit;

for i=1:size(FIMR.fimResultsOne,1)
    fimPrior = diag(1./covPrior0);
    
    FIM = NCells(1)*FIMR.fimResultsOne{i,iPDO}{2};
    FIM = diag(FIMR.parsets{iPDO}(i,:))*FIM*diag(FIMR.parsets{iPDO}(i,:));
    FIM = FIM*log(10)^2+fimPrior;
    detFIM_Total_Best_0(i,:) = det(FIM(pars,pars)^-1);
    sigsFIM_Total_Best_0(i,:) = diag(FIM(pars,pars)^(-1));
    

    FIM = NCells(1)*FIMR.fimResultsOne{i,iPDO}{2}+...
        NCells(2)*FIMR.fimResultsOne{i,iPDO}{6};
    FIM = diag(FIMR.parsets{iPDO}(i,:))*FIM*diag(FIMR.parsets{iPDO}(i,:));
    FIM = FIM*log(10)^2+fimPrior;
    detFIM_Total_Best_0_18(i,:) = det(FIM(pars,pars)^-1);
    sigsFIM_Total_Best_0_18(i,:) = diag(FIM(pars,pars)^(-1));

    FIM = NCells(1)*FIMR.fimResultsOne{i,iPDO}{2}+...
        NCells(3)*FIMR.fimResultsOne{i,iPDO}{53};
    FIM = diag(FIMR.parsets{iPDO}(i,:))*FIM*diag(FIMR.parsets{iPDO}(i,:));
    FIM = FIM*log(10)^2+fimPrior;
    detFIM_Total_Best_0_300(i,:) = det(FIM(pars,pars)^-1);
    sigsFIM_Total_Best_0_300(i,:) = diag(FIM(pars,pars)^(-1));
    if i==1
        FIM_0_300=FIM;
    end

    FIM = NCells(1)*FIMR.fimResultsOne{i,iPDO}{2}+...
        NCells(2)*FIMR.fimResultsOne{i,iPDO}{6}+...
        NCells(3)*FIMR.fimResultsOne{i,iPDO}{53};
    FIM = diag(FIMR.parsets{iPDO}(i,:))*FIM*diag(FIMR.parsets{iPDO}(i,:));
    FIM = FIM*log(10)^2+fimPrior;
    detFIM_Total_Best_0_18_300(i,:) = det(FIM(pars,pars)^-1);
    sigsFIM_Total_Best_0_18_300(i,:) = diag(FIM(pars,pars)^(-1));
end

% Final_Total = [detFIM_Total_Best_0_18(1,:),...
%     detFIM_Total_Best_0_300(1,:),detFIM_Total_Best_0_18_300(1,:)];
Final_Total = [mean((detFIM_Total_Best_0_18),1),...
    mean((detFIM_Total_Best_0_300),1),mean((detFIM_Total_Best_0_18_300),1)];
Sigs_Total = [mean((sigsFIM_Total_Best_0_18),1),...
    mean((sigsFIM_Total_Best_0_300),1),mean((sigsFIM_Total_Best_0_18_300),1)];
Final_STD = [std((detFIM_Total_Best_0_18),[],1),...
    std((detFIM_Total_Best_0_300),[],1),std((detFIM_Total_Best_0_18_300),[],1)];
end

function [BIC] = makeFiguresForPDOs(vars,FN,NCells,kPars)

PDOStats = FittingFunctionsFISHTrue(1,[FN,'_0_300'],vars);

figure(1);
for i=1:3
    subplot(3,1,i)
    BIC{i} = [sum(PDOStats.loglPred([1:3]+(i-1)*3,[1,3]),2),PDOStats.loglPred([1:3]+(i-1)*3,2)];
    BIC{i}(:,1) = -2*BIC{i}(:,1)+log(NCells(1)+NCells(3))*kPars;
    BIC{i}(:,2) = -2*BIC{i}(:,2)+log(NCells(2))*[3;4;7];
end

%% Make PDO Figures
vars.doFit = 0;
vars.makePlots = true;
FittingFunctionsFISHTrue(1,[FN,'_0_300'],vars);

%%         Reformat for manuscript
legs = {'FISH Spots','MCP-GFP Spots','PDO \times FISH';...
    'FISH Spots','MCP-GFP Spots','PDO \times FISH';...
    'FISH Spots','MCP-GFP Spots','PDO \times FISH';...
    'FISH Spots','FISH Intens','PDO \times FISH';...
    'FISH Spots','FISH Intens','PDO \times FISH';...
    'FISH Spots','FISH Intens','PDO \times FISH';...
    'FISH Spots','MCP-GFP Intens','PDO \times FISH';...
    'FISH Spots','MCP-GFP Intens','PDO \times FISH';...
    'FISH Spots','MCP-GFP Intens','PDO \times FISH'};

for icontourFig = [3,5,9]
    f = figure(icontourFig);
    fa = gca;
    set(f,'Position',[200*icontourFig   278   338   269])
    set(gca,'xlim',[0,500])
    set(gca,'XLabel',[],'YLabel',[])
    ylim = get(gca,'ylim');

    f2 = figure(icontourFig+10);clf; set(f2,'Position',[200*icontourFig   278   338   269])
    f2a = gca;
    copyobj(fa.Children([3,4]),f2a);
    set(gca,'xlim',[0,500],'xlim',[0,500],'ylim',ylim,'FontSize',15)
    colorbar

    f2 = figure(icontourFig+20);clf; set(f2,'Position',[200*icontourFig   278   338   269])
    f2a = gca;
    copyobj(fa.Children([2,4]),f2a);
    set(gca,'xlim',[0,500],'ylim',ylim,'FontSize',15)
    colorbar

    f2 = figure(icontourFig+30);clf; set(f2,'Position',[200*icontourFig   278   338   269])
    f2a = gca;
    copyobj(fa.Children([1,4]),f2a);
    set(gca,'xlim',[0,500],'ylim',ylim,'FontSize',15)
    colorbar


    for iTime = 1:3
        f = figure(100*iTime+icontourFig); set(f,'Position',[200*icontourFig  678   338   269])
        set(gca,'XLabel',[],'YLabel',[])
        legend(legs(icontourFig,:),'Location','southeast')
        if icontourFig==3
            set(gca,'xlim',[0,500],'ylim',[0,1]);
        elseif icontourFig==5
            set(gca,'xlim',[0,2500]);
        elseif icontourFig==9
            set(gca,'xlim',[0,3000],'ylim',[0,1]);
        end

        f = figure(1000+100*iTime+icontourFig); set(f,'Position',[200*icontourFig    678   338   269])
        set(gca,'XLabel',[],'YLabel',[])
        legend(legs(icontourFig,:),'Location','southeast')
        if icontourFig==3
            set(gca,'xlim',[0,500],'ylim',[0,0.012]);
        elseif icontourFig==5
            set(gca,'xlim',[0,2500],'ylim',[0,0.012]);
        elseif icontourFig==9
            set(gca,'xlim',[0,3000],'ylim',[0,0.012]);
        end
    end
end
end

function [Fitresults] = makeFitFigures(vars,FNext)

Fitresults = FittingFunctionsFISHTrue(32,FNext,vars);

%%         Reformat for Manuscript
clear hOrig1 hOrig2
for i = [1:5]
    for j = [1,2,3]
        f = figure(i);
        set(f,'position',[560   694   811   253])
        subplot(1,3,j);
        hOrig1{i,j} = gca;
        switch i
            case {1,2,3}
                set(gca,'xlim',[-5,600],'ylim',[0,1.05])
            case {4,5}
                set(gca,'xlim',[-5,3000],'ylim',[0,1.05])
        end
        
        f = figure(i+10);
        set(f,'position',[560   294   811   253])
        subplot(1,3,j);
        hOrig2{i,j} = gca;
        switch i
            case {1,2,3}
                set(gca,'xlim',[-5,600],'ylim',[0,0.2])
            case {4,5}
                set(gca,'xlim',[-5,3000],'ylim',[0,0.2])
        end
    end
end

%
iarr = [1,3,4,5];
jarr = [1,3,2];
for k=1:2
    switch k
        case 1
            nf = figure(117);clf;
            hOrig = hOrig1;
        case 2
            nf = figure(118);clf;
            hOrig = hOrig2;
    end

    for i=1:4
        for j = 1:3
            figure(nf)
            subplot(4,4,(i-1)*4+j)
            na = gca;
            if i==1&&j<3
                copyobj(hOrig{iarr(i),jarr(j)}.Children(1:2),na)
            elseif j<3
                copyobj(hOrig{iarr(i),jarr(j)}.Children(3:4),na)
            elseif i==1&&j==3
                subplot(4,4,(i-1)*4+4)
                na = gca;
                copyobj(hOrig{iarr(i),jarr(3)}.Children(1:2),na)
                na.Children(1).Color=[0.494 0.184 0.556];
                na.Children(2).Color=[ 0.929 0.694 0.125];
            else
                subplot(4,4,(i-1)*4+3)
                na = gca;
                copyobj(hOrig{iarr(i),jarr(3)}.Children(3:4),na)
                subplot(4,4,(i-1)*4+4)
                na.Children(1).Color=[0.494 0.184 0.556];
                na.Children(2).Color=[ 0.929 0.694 0.125];
                na = gca;
                copyobj(hOrig{iarr(i),jarr(3)}.Children(1:2),na)
                na.Children(1).Color=[0.494 0.184 0.556];
                na.Children(2).Color=[ 0.929 0.694 0.125];
            end
        end
    end
    for i=1:4
        for j = 1:4
            subplot(4,4,(i-1)*4+j)
            na = gca;
            if i<=2||j==4
                set(na,'fontsize',15,'xlim',[0,500],'xtick',[0,250,500])
                if k==2&&j==2
                    set(gca,'YLim',[0,1])
                elseif k==2&&j==1
                    set(gca,'YLim',[0,0.23])
                end                    
            else
                set(na,'fontsize',15,'xlim',[0,2500])
                if k==2
                    set(gca,'YLim',[0,0.1])
                end

            end
        end
    end
end
end

function makeMHandFimPlots(vars,FN,FNext,adjustLims,ax1,ax4)
arguments
    vars
    FN
    FNext
    adjustLims=true;
    ax1 = [-4.5,-3.5; -3,1; -1,3; -2.5,-1.5];
    ax4 = [-6.5,-3.5; -3,1; -1,3; -3,-0];
end
    
vars.mleScatterPlots = true;
vars.covWithPrior=true;
FittingFunctionsFISHTrue(4,FNext,vars)

load([FN,'_0_18_300'],'ModelZero');
x = log10([ModelZero{1}.parameters{2:end,2}]);
npar = length(x);
for k = 2:3:14
    figure(k)
    for i=1:npar-1
        for j = i:npar-1
            subplot(npar-1,npar-1,(i-1)*(npar-1)+j)
            plot([-10,10],x(i)*[1,1],'k--')
            plot(x(j+1)*[1,1],[-10,10],'k--')
        end
    end
end

if ~adjustLims
    return
end

%%         Adjust figures to zoom in on support and denote correct paramters.
for j=[2,5,8,11,14]
    if j<14
        ax=ax1;
    else
        ax=ax4;
    end
figure(j)
m = size(ax,1);
for i1=1:m-1
    for i2 = i1:m-1
        subplot(m-1,m-1,(i1-1)*(m-1)+i2);
        set(gca,'ylim',ax(i1,:))
        set(gca,'xlim',ax(i2+1,:))
    end
end
%         
% for i=1:npar-1
%     subplot(npar-1,npar-1,i)
%     set(gca,'ylim',ax(1,:))
% end
% subplot(npar-1,npar-1,1)
% set(gca,'xlim',ax(2,:))
% for i=npar+1:2*(npar-1)
%     subplot(npar-1,npar-1,i)
%     set(gca,'ylim',ax(2,:))
% end
% if npar>=4
%     for i=npar-1:npar-1:(npar-1)^2
%         subplot(npar-1,npar-1,i)
%         set(gca,'xlim',ax(4,:))
%     end
%     subplot(npar-1,npar-1,9)
%     set(gca,'ylim',ax(3,:))
% end
% for i=2:npar-1:npar+1
%     subplot(npar-1,npar-1,i)
%     set(gca,'xlim',ax(3,:))
% end
end
end

function makeTexParTable(FNext)
load(FNext);
Atxt = importdata('InitParTableTex2.txt');
va = ['O';'W';'B';'G'];
vn = [2:5];
pdo = [1,3,4,5];

for i=1:4
    for j=1:4
        Atxt = strrep(Atxt,[va(i),'1',num2str(j)],num2str(log10(ModelZero{pdo(j)}.parameters{vn(i),2}),'%0.3f'));

        smplDone = chainResults{pdo(j)}.mhSamples(:,:);
        valDone = chainResults{pdo(j)}.mhValue;
        smplDone = smplDone(valDone~=0,:);
        covMH = cov(smplDone/log(10));
        Atxt = strrep(Atxt,[va(i),'2',num2str(j)],num2str(sqrt(covMH(i,i)),'%0.3f'));

        covFIM = (FIMZeroLog{pdo(j)}+diag([2 1 1 1]*log(10))^(-1))^(-1);
        covFIM = covFIM/(log(10)^2);

%         covFIM = covLog{pdo(j)}/(log(10)^2);
        Atxt = strrep(Atxt,[va(i),'3',num2str(j)],num2str(sqrt(covFIM(i,i)),'%0.3f'));

    end
end
clc
for i=1:length(Atxt)
    disp(Atxt{i})
end
end

function makeFimExptDesign(vars,FN,lims)
FittingFunctionsFISHTrue(6,[FN,'_0_300'],vars)
% set(gca,'xlim',[0,800])
f3 =gcf; set(f3,'Name','FIM vs t2, Fit on t=0,300')
a3=gca;

FittingFunctionsFISHTrue(6,[FN,'_0_18_300'],vars) 
% set(gca,'xlim',[0,800])
f4 =gcf; set(f4,'Name','FIM vs t2, Fit on t=0,18,300')
a4=gca;
%%         Reformat for manuscript
f = figure(1112); clf;
set(f,'Position',[ 855   524   889   196]);
for j=1:4
    n{j} = subplot(1,4,5-j);
    copyobj(a3.Children(j+0:4:16),n{j})
    set(n{j}.Children(4),'FaceColor','m'); 
    set(n{j}.Children(3),'Color','k');
    set(n{j}.Children(2),'Color','m'); 
    set(n{j}.Children(1),'Color','k'); n{j}.Children(1).Color(4) = 0.25;
    
    copyobj(a4.Children(j+0:4:16),n{j})
    set(n{j}.Children(4),'FaceColor','c','FaceAlpha',0.5); 
    set(n{j}.Children(3),'Color','k')
    set(n{j}.Children(2),'Color','k'); n{j}.Children(2).Color(4) = 0.25;
    set(n{j}.Children(1),'Color','c'); 

    set(n{j},'FontSize',15)
%     legend(n{j}.Children([8,4]),{'0,300','0,18,300'},'Location','northeast','FontSize',12)
end

% set(n{4},'xlim',[0,3.2],'ylim',[lims(1,1),lims(1,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'})
% set(n{3},'xlim',[0,3.2],'ylim',[lims(2,1),lims(2,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'})
% set(n{2},'xlim',[0,3.2],'ylim',[lims(3,1),lims(3,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'})
% set(n{1},'xlim',[0,3.2],'ylim',[lims(4,1),lims(4,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'})
for i=1:4
    for j = 1:(lims(i,2)-lims(i,1)+1)
        limtext{i}{j} = ['10^{',num2str(lims(i,1)+j-1),'}'];
    end
end

set(n{4},'xlim',[0,3.2],'ylim',[lims(1,1),lims(1,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'},...
    'ytick',[lims(1,1):lims(1,2)],'yticklabel',limtext{1}(:)')
set(n{3},'xlim',[0,3.2],'ylim',[lims(2,1),lims(2,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'},...
    'ytick',[lims(2,1):lims(2,2)],'yticklabel',limtext{2}(:)')
set(n{2},'xlim',[0,3.2],'ylim',[lims(3,1),lims(3,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'},...
    'ytick',[lims(3,1):lims(3,2)],'yticklabel',limtext{3}(:)')
set(n{1},'xlim',[0,3.2],'ylim',[lims(4,1),lims(4,2)],'xtick',[0,1,2,3],'xticklabel',{'10^0','10^1','10^2','10^3'},...
    'ytick',[lims(4,1):lims(4,2)],'yticklabel',limtext{4}(:)')

end

function [detCovFIM,detCovMH,FIM_0_300] = makeFimExpComparisons(FN,NCells,vars)
% FIM using the final parameters
FIMR_Final = load([FN,'_0_18_300'],'fimResultsOne','parsets');

[Final_FISH,FinalsigsFISH,Final_FISH_sd] = evaluateFIMDesigns(FIMR_Final,1,NCells,vars);
[Final_MCPspots,FinalsigsMCPspots,Final_MCPspots_sd] = evaluateFIMDesigns(FIMR_Final,3,NCells,vars);
[Final_smiFISHIntens,FinalsigsFISHIntens,Final_smiFISHIntens_sd] = evaluateFIMDesigns(FIMR_Final,4,NCells,vars);
[Final_smiMCPIntens,FinalsigsMCPIntens,Final_smiMCPIntens_sd] = evaluateFIMDesigns(FIMR_Final,5,NCells,vars);

% FIM using the initial parameters (after fitting t=0,300)
FIMR_w_0_300 = load([FN,'_0_300'],'fimResultsOne','parsets');
[w_0_300_FISH,w_0_300_sigsFISH,w_0_300_FISH_sd,FIM_0_300{1}] = evaluateFIMDesigns(FIMR_w_0_300,1,NCells,vars);
[w_0_300_MCPspots,w_0_300_sigsMCPspots,w_0_300_MCPspots_sd,FIM_0_300{2}] = evaluateFIMDesigns(FIMR_w_0_300,3,NCells,vars);
[w_0_300_smiFISHIntens,w_0_300_sigsFISHIntens,w_0_300_smiFISHIntens_sd,FIM_0_300{3}] = evaluateFIMDesigns(FIMR_w_0_300,4,NCells,vars);
[w_0_300_smiMCPIntens,w_0_300_sigsMCPIntens,w_0_300_smiMCPIntens_sd,FIM_0_300{4}] = evaluateFIMDesigns(FIMR_w_0_300,5,NCells,vars);

% FIM using the initial parameters (after fitting t=0,18)
FIMR_w_0_18 = load([FN,'_0_18'],'fimResultsOne','parsets');
[w_0_18_FISH,w_0_18_sigsFISH] = evaluateFIMDesigns(FIMR_w_0_18,1,NCells,vars);
[w_0_18_MCPspots,w_0_18_sigsMCPspots] = evaluateFIMDesigns(FIMR_w_0_18,3,NCells,vars);
[w_0_18_smiFISHIntens,w_0_18_sigsFISHIntens] = evaluateFIMDesigns(FIMR_w_0_18,4,NCells,vars);
[w_0_18_smiMCPIntens,w_0_18_sigsMCPIntens] = evaluateFIMDesigns(FIMR_w_0_18,5,NCells,vars);

figure(1); clf
set(gcf,'position',[560   565   889   382])
subplot(2,4,1)
bar([1,2,3]-.23,[w_0_300_FISH]',.44,'m'); hold on;
bar([1,2,3]+.23,[Final_FISH]',.44,'c')
errorbar([1,2,3]+.23,Final_FISH',Final_FISH_sd','c.','linewidth',3)
errorbar([1,2,3]-.23,w_0_300_FISH',w_0_300_FISH_sd','m.','linewidth',3)

subplot(2,4,2)
bar([1,2,3]-.23,[w_0_300_MCPspots]',.44,'m'); hold on;
bar([1,2,3]+.23,[Final_MCPspots]',.44,'c')
errorbar([1,2,3]+.23,Final_MCPspots',Final_MCPspots_sd','c.','linewidth',3)
errorbar([1,2,3]-.23,w_0_300_MCPspots',w_0_300_MCPspots_sd','m.','linewidth',3)

subplot(2,4,3)
bar([1,2,3]-.23,[w_0_300_smiFISHIntens]',.44,'m'); hold on;
bar([1,2,3]+.23,[Final_smiFISHIntens]',.44,'c')
errorbar([1,2,3]+.23,Final_smiFISHIntens',Final_smiFISHIntens_sd','c.','linewidth',3)
errorbar([1,2,3]-.23,w_0_300_smiFISHIntens',w_0_300_smiFISHIntens_sd','m.','linewidth',3)

subplot(2,4,4)
bar([1,2,3]-.23,[w_0_300_smiMCPIntens]',.44,'m'); hold on;
bar([1,2,3]+.23,[Final_smiMCPIntens]',.44,'c') 
errorbar([1,2,3]+.23,Final_smiMCPIntens',Final_smiMCPIntens_sd','c.','linewidth',3)
errorbar([1,2,3]-.23,w_0_300_smiMCPIntens',w_0_300_smiMCPIntens_sd','m.','linewidth',3)

detCovFIM = [w_0_300_FISH(2),w_0_300_MCPspots(2),w_0_300_smiFISHIntens(2),w_0_300_smiMCPIntens(2)];

%         Add Spread of MH results before and After 2nd Experiment
load([FN,'_0_18'],'chainResults');
for i=[1,3,4,5]
    detCovarTotal(i,1) = det(cov(chainResults{i}.mhSamples/log(10)));
    sigsTotal(i,1:size(vars.modelVarsToFit,2)) = diag(cov(chainResults{i}.mhSamples/log(10)));
end

load([FN,'_0_300'],'chainResults');
for i=[1,3,4,5]
    detCovarTotal(i,2) = det(cov(chainResults{i}.mhSamples/log(10)));
    sigsTotal(i,size(vars.modelVarsToFit,2)+1:size(vars.modelVarsToFit,2)*2) = diag(cov(chainResults{i}.mhSamples/log(10)));
end
load([FN,'_0_18_300'],'chainResults');
for i=[1,3,4,5]
    detCovarTotal(i,3) = det(cov(chainResults{i}.mhSamples/log(10)));
    sigsTotal(i,size(vars.modelVarsToFit,2)*2+1:size(vars.modelVarsToFit,2)*3) = diag(cov(chainResults{i}.mhSamples/log(10)));
end

detCovMH = detCovarTotal([1,3:5],2)';
% detCovMH = detCovarTotal(:,2)';

figure(1)
subplot(2,4,5)
bar(detCovarTotal(1,:),'b')

subplot(2,4,6)
bar(detCovarTotal(3,:),'b')

subplot(2,4,7)
bar(detCovarTotal(4,:),'b')

subplot(2,4,8)
bar(detCovarTotal(5,:),'b')

for i=1:4
    subplot(2,4,i)
    set(gca,'yscale','log','ylim',[1e-13,3e-4],'xlim',[0 4],'xtick',[],'ytick',10.^[-12:2:-4],'FontSize',15)
    
    subplot(2,4,i+4)
    set(gca,'yscale','log','ylim',[1e-13,3e-4],'xlim',[0 4],'xticklabel',{'0,18','0,300','0,18,300'},'ytick',10.^[-12:2:-4],'FontSize',15)
end
%%   Make scatter plot of expected and actual COV determinants.
f = figure(22); clf;
set(f,'position',[ 1465         988         319         236])
AA = [detCovarTotal(1,:),detCovarTotal(3,:),detCovarTotal(4,:),detCovarTotal(5,:)];
BB_0_300 = [[w_0_300_FISH],[w_0_300_MCPspots],[w_0_300_smiFISHIntens],[w_0_300_smiMCPIntens]];
BB_0_18 = [[w_0_18_FISH],[w_0_18_MCPspots],[w_0_18_smiFISHIntens],[w_0_18_smiMCPIntens]];
BB_Final = [[Final_FISH],[Final_MCPspots],[Final_smiFISHIntens],[Final_smiMCPIntens]];

% BB2 = [[w_0_18_FISH(1),w_0_300_FISH(2),Final_FISH(3)],...
%     [w_0_18_MCPspots(1),w_0_300_MCPspots(2),Final_MCPspots(3)],...
%     [w_0_18_smiFISHIntens(1),w_0_300_smiFISHIntens(2),Final_smiFISHIntens(3)],...
%     [w_0_18_smiMCPIntens(1),w_0_300_smiMCPIntens(2),Final_smiMCPIntens(3)]];
scatter(log10(BB_0_18),log10(AA),200,'m^','filled'); hold on
scatter(log10(BB_0_300),log10(AA),200,'bo','filled'); hold on
scatter(log10(BB_Final),log10(AA),200,'ks','filled'); hold on
plot([-12 -3],[-12,-3],'k--')
set(gca,'fontsize',15,'xlim',[-12 -4],'ylim',[-12 -4])
legend({'Fit t \in \{0,18\}','Fit t \in \{0,300\}','Fit t \in \{0,18,300\}'},'Location','Northwest')

%%   Make scatter plots of expected and actual parameter variances
figure(23); clf;
set(gcf,'Position',[561   747   889   160]);
mkr = ['o';'s';'d';'v']; % Markers correspond to different parameters.
col = ['k','m','c','b']; % Colors correspond to different experiment designs.

npar = size(vars.modelVarsToFit,2);
for i=1:npar
    %     for j=1:3
    subplot(1,5,1)

    plot(w_0_18_sigsFISH(i),sigsTotal(1,i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(1)); hold on
    plot(w_0_300_sigsFISH(npar+i),sigsTotal(1,npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(2)); hold on
    plot(FinalsigsFISH(2*npar+i),sigsTotal(1,2*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(3)); hold on

    subplot(1,5,2)
    plot(w_0_18_sigsMCPspots(i),sigsTotal(3,i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(1)); hold on
    plot(w_0_300_sigsMCPspots(npar+i),sigsTotal(3,npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(2)); hold on
    plot(FinalsigsMCPspots(2*npar+i),sigsTotal(3,2*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(3)); hold on
    %
    subplot(1,5,3)
    plot(w_0_18_sigsFISHIntens(i),sigsTotal(4,i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(1)); hold on
    plot(w_0_300_sigsFISHIntens(npar+i),sigsTotal(4,npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(2)); hold on
    plot(FinalsigsFISHIntens(2*npar+i),sigsTotal(4,2*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(3)); hold on

    subplot(1,5,4)
    plot(w_0_18_sigsMCPIntens(i),sigsTotal(5,i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(1)); hold on
    plot(w_0_300_sigsMCPIntens(npar+i),sigsTotal(5,npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(2)); hold on
    plot(FinalsigsMCPIntens(2*npar+i),sigsTotal(5,2*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(3)); hold on

        subplot(1,5,5)
        plot(w_0_18_sigsFISH((1-1)*npar+i),sigsTotal(1,(1-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(1)); hold on
        plot(w_0_300_sigsFISH((2-1)*npar+i),sigsTotal(1,(2-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(1)); hold on
        plot(FinalsigsFISH((3-1)*npar+i),sigsTotal(1,(3-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(1)); hold on

        plot(w_0_18_sigsMCPspots((1-1)*npar+i),sigsTotal(3,(1-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(2)); hold on
        plot(w_0_300_sigsMCPspots((2-1)*npar+i),sigsTotal(3,(2-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(2)); hold on
        plot(FinalsigsMCPspots((3-1)*npar+i),sigsTotal(3,(3-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(2)); hold on

        plot(w_0_18_sigsFISHIntens((1-1)*npar+i),sigsTotal(4,(1-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(3)); hold on
        plot(w_0_300_sigsFISHIntens((2-1)*npar+i),sigsTotal(4,(2-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(3)); hold on
        plot(FinalsigsFISHIntens((3-1)*npar+i),sigsTotal(4,(3-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(3)); hold on

        plot(w_0_18_sigsMCPIntens((1-1)*npar+i),sigsTotal(5,(1-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(4)); hold on
        plot(w_0_300_sigsMCPIntens((2-1)*npar+i),sigsTotal(5,(2-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(4)); hold on
        plot(FinalsigsMCPIntens((3-1)*npar+i),sigsTotal(5,(3-1)*npar+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(4)); hold on

end

for i=1:5
    subplot(1,5,i)
    set(gca,'YScale','log','xscale','log','fontsize', 15)
    set(gca,'xlim',10.^[-3.5 1],'ylim',10.^[-3.5 1],'xtick',10.^[-2:2:1],'ytick',10.^[-2:2:1])
    set(gca,'XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off')
    plot([1e-4 1e2],[1e-4 1e2],'k--')
end
end