clear all
addpath('../')
TMP = SSIT(); % Load definition of the SSIT class for later use.
clear TMP
%% Make PDO Figures
vars.doFit = 0;
vars.modelVarsToFit = [1,2,3,4];
vars.pdoTimes = [0,300];
FN = 'FinalFit';
FittingFunctionsCoLocalized(1,[FN],vars)

%%     Reformat for manuscript
f = figure(2); set(f,'Position',[560   678   338   269])
legend('PDO','0 min','18 min','300 min','Location','northwest')
set(gca,'XLabel',[],'YLabel',[])

f = figure(22); set(f,'Position',[560   678   338   269])
legend('PDO','0 min','18 min','300 min','Location','northwest')
set(gca,'XLabel',[],'YLabel',[])

f = figure(105); set(f,'Position',[560   678   338   269])
set(gca,'XLabel',[],'YLabel',[])
legend('Location','southeast')
f = figure(125); set(f,'Position',[560   678   338   269])
set(gca,'XLabel',[],'YLabel',[])
legend('Location','southeast')
f = figure(106); set(f,'Position',[560   678   338   269])
set(gca,'XLabel',[],'YLabel',[],'ylim',[0,0.01])
legend('Location','northeast')
f = figure(126); set(f,'Position',[560   678   338   269])
set(gca,'XLabel',[],'YLabel',[],'ylim',[0,0.01])
legend('Location','northeast')

for i=2:3
    f = figure(100*i+5); set(f,'Position',[560   678   338   269])
    legend('Location','southeast')
    f = figure(100*i+25); set(f,'Position',[560   678   338   269])
    legend('Location','southeast')
    f = figure(100*i+6); set(f,'Position',[560   678   338   269])
    legend('Location','northeast')
    f = figure(100*i+26); set(f,'Position',[560   678   338   269])
    legend('Location','northeast')
end

%% Make Initial Fit Figures (fitting 0,300 only)
close all force
vars.doFit = 0;
vars.modelVarsToFit = [2:5];
FittingFunctionsCoLocalized(32,[FN,'_0_300'],vars)
close(4);close(3)

clear h1 h2 h20 h10
for i = 1:2
    if i==1;k=1;else;k=3;end
    figure(2);
    subplot(1,3,k);
    h2{i}=gca; 
    
    figure(20); 
    subplot(1,3,k);
    h20{i}=gca; 

    figure(1);
    subplot(1,3,k);
    h1{i}=gca; 
    
    figure(10); 
    subplot(1,3,k);
    h10{i}=gca; 
end
%%     Reformat for manuscript
f = figure(1111); clf; set(f,'position',[153   422   690   710])
for i=1:12
    n{i} = subplot(3,4,i);
end
for i=1:2
    copyobj(h10{i}.Children,n{i});
    copyobj(h1{i}.Children,n{i+2});
    copyobj(h20{i}.Children,n{4+i});
    copyobj(h2{i}.Children,n{6+i});
    copyobj(h20{i}.Children,n{8+i});
    copyobj(h2{i}.Children,n{10+i});
end

for i=[2:4,6:8,10:12]
    set(n{i},'YTickLabel',[])
end

for i=1:4
    set(n{i},'FontSize',14,'xlim',[0,1000],'YLim',[0,1.01])
end
for i=5:12
    set(n{i},'FontSize',14,'xlim',[0,1000],'YLim',[0,0.1])
end

legend(n{1}.Children([8,6,7,5,3,1]),{'Total Data','FISH Data',...
    'Total Fit','FISH Fit','PDO Fit',...
    'PDO Prediction'},'Location','southeast','FontSize',9)
legend(n{3}.Children([8,6,7,5,3,1]),{'Total Data','MCP Data',...
    'Total Fit','MCP Fit','PDO Fit',...
    'PDO Prediction'},'Location','southeast','FontSize',9)

cols([1,3,7,5,8,6]) = ['ccmmkk'];
for tim = 1:2
    for gene=1:2
        for msmt = 1:2
            switch gene
                case 2
                    nm = 'MCP';
                case 1
                    nm = 'FISH';
            end
            switch msmt
                case 1
                    lineShow = [8,7,1];
                case 2
                    lineShow = [6,5,3];
            end
            k = 4+(2-msmt)*4+(gene-1)*2+tim;
            for i=1:8
                n{k}.Children(i).Visible = 'off';
            end
            for i=lineShow
                n{k}.Children(i).Visible = 'on';
                set(n{k}.Children(i),'color',cols(i))
                set(n{k}.Children(i),'LineStyle','-')
            end
            set(n{k},'ylim',[0,0.1],'xlim',[0,1000],'fontsize',14)
            if msmt==2&&tim==1
                legend(n{k}.Children(lineShow),{[nm,' Data'],[nm,' Fit'],['PDO Fit']},'Location','northeast','FontSize',9)
            elseif tim==1
                legend(n{k}.Children(lineShow),{'Total Data','Total Fit',['PDO Prediction']},'Location','northeast','FontSize',9)
            end
            if tim==2||gene==2
                set(n{k},'YTickLabel',[])
            end
        end
    end
end

%% Make Initial MH Scatter Plots (fitting 0,300 only)
close all
vars.doFit = 0;
vars.priorScale = 1;
vars.mleScatterPlots = true;
vars.modelVarsToFit = [2:5];
FittingFunctionsCoLocalized(4,[FN,'_0_300'],vars)

%%         Adjust figures to zoom in on support
ax = [-4,-3.4; -2,-0.8; 1.0,2.2; -2.4,-2.1];
for j=[2,11,14]
figure(j)
for i=1:3
    subplot(3,3,i)
    set(gca,'ylim',ax(1,:))
end
subplot(3,3,1)
set(gca,'xlim',ax(2,:))
for i=5:6
    subplot(3,3,i)
    set(gca,'ylim',ax(2,:))
end
for i=3:3:9
    subplot(3,3,i)
    set(gca,'xlim',ax(4,:))
end
for i=2:3:5
    subplot(3,3,i)
    set(gca,'xlim',ax(3,:))
end
subplot(3,3,9)
set(gca,'ylim',ax(3,:))
end

%% Make Tex Table for Initial parameters;
load([FN,'_0_300']);
Atxt = importdata('InitParTableTex.txt');
va = ['O';'W';'B';'G'];
vn = [2:5];
pdo = [1,4,5];

for i=1:4
    for j=1:3
        Atxt = strrep(Atxt,[va(i),'1',num2str(j)],num2str(log10(ModelZero{pdo(j)}.parameters{vn(i),2}),'%0.3f'));

        smplDone = chainResults{pdo(j)}.mhSamples(:,:);
        valDone = chainResults{pdo(j)}.mhValue;
        smplDone = smplDone(valDone~=0,:);
        covMH = cov(smplDone/log(10));
        Atxt = strrep(Atxt,[va(i),'2',num2str(j)],num2str(sqrt(covMH(i,i)),'%0.3f'));

        covFIM = covLog{pdo(j)}/(log(10)^2);
        Atxt = strrep(Atxt,[va(i),'3',num2str(j)],num2str(sqrt(covFIM(i,i)),'%0.3f'));

    end
end
Atxt

%% Make Plot of FIM vs time of third experiment
close all
vars.doFit = 1;
vars.priorScale = 1;
vars.FIMShowErrorBars = true;
vars.FIMnCells = [155,0,66];

FittingFunctionsCoLocalized(6,[FN,'_0_300'],vars)
set(gca,'xlim',[0,800])
f3 =gcf; set(f3,'Name','FIM vs t2, Fit on t=0,300')
a3=gca;

FittingFunctionsCoLocalized(6,[FN,'_0_18_300'],vars) 
set(gca,'xlim',[0,800])
f4 =gcf; set(f4,'Name','FIM vs t2, Fit on t=0,18,300')
a4=gca;
%%       Reformat for manuscript
f = figure(1112); clf;
set(f,'Position',[382   250   771   233]);
for j=1:3
    n{j} = subplot(1,3,4-j);
    copyobj(a3.Children(j+0:3:12),n{j})
    set(n{j}.Children(4),'Color','m')
    set(n{j}.Children(3),'Color','k')
    set(n{j}.Children(2),'Color','m')
    set(n{j}.Children(1),'Color','k')
    
    copyobj(a4.Children(j+0:3:12),n{j})
    set(n{j}.Children(4),'Color','c')
    set(n{j}.Children(3),'Color','k')
    set(n{j}.Children(2),'Color','k')
    set(n{j}.Children(1),'Color','c')

    set(n{j},'FontSize',15,'xlim',[0,800],'ylim',[-12,-7.5])
    legend(n{j}.Children([8,4]),{'0,300','0,18,300'},'Location','northeast','FontSize',12)
end

set(n{1},'ylim',[-10.5,-8])
set(n{2},'ylim',[-10.5,-8])
set(n{3},'ylim',[-10.5,-8])


%% Compare FIM for different experiment designs.
close all

% FIM using the final parameters
FIMR_Final = load([FN,'_0_18_300'],'fimResultsOne','parsets');
NCells = [155 106 66];

[Final_Total,FinalsigsTotal] = evaluateFIMDesigns(FIMR_Final,1,NCells,[2:5]);
[Final_MCP,FinalsigsMCP] = evaluateFIMDesigns(FIMR_Final,4,NCells,[2:5]);
[Final_smiFISH,FinalsigsFISH] = evaluateFIMDesigns(FIMR_Final,5,NCells,[2:5]);

% FIM using the initial parameters (after fitting t=0,300)
FIMR_Initial = load([FN,'_0_300'],'fimResultsOne','parsets');
figure(1); clf
[Initial_Total,InitialsigsTotal] = evaluateFIMDesigns(FIMR_Initial,1,NCells,[2:5]);
subplot(2,3,1)
bar([1,2,3]-.23,1./[Initial_Total]',.46,'m'); hold on;
bar([1,2,3]+.23,1./[Final_Total]',.46,'c')

[Initial_MCP,InitialsigsMCP] = evaluateFIMDesigns(FIMR_Initial,4,NCells,[2:5]);
subplot(2,3,2)
bar([1,2,3]-.23,1./[Initial_MCP]',.46,'m'); hold on;
bar([1,2,3]+.23,1./[Final_MCP]',.46,'c')

[Initial_smiFISH,InitialsigsFISH] = evaluateFIMDesigns(FIMR_Initial,5,NCells,[2:5]);
subplot(2,3,3)
bar([1,2,3]-.23,1./[Initial_smiFISH]',.46,'m'); hold on;
bar([1,2,3]+.23,1./[Final_smiFISH]',.46,'c')
%%         Add Spread of MH results before and After 2nd Experiment
load([FN,'_0_18'],'chainResults');
for i=[1,4,5]
    detCovarTotal(i,1) = det(cov(chainResults{i}.mhSamples));
    sigsTotal(i,1:4) = diag(cov(chainResults{i}.mhSamples));
end

load([FN,'_0_300'],'chainResults');
for i=[1,4,5]
    detCovarTotal(i,2) = det(cov(chainResults{i}.mhSamples));
    sigsTotal(i,5:8) = diag(cov(chainResults{i}.mhSamples));
end
load([FN,'_0_18_300'],'chainResults');
for i=[1,4,5]
    detCovarTotal(i,3) = det(cov(chainResults{i}.mhSamples));
    sigsTotal(i,9:12) = diag(cov(chainResults{i}.mhSamples));
end

figure(1)
subplot(2,3,4)
bar(detCovarTotal(1,:),'b')

subplot(2,3,5)
bar(detCovarTotal(4,:),'b')

subplot(2,3,6)
bar(detCovarTotal(5,:),'b')

for i=1:3
    subplot(2,3,i)
    set(gca,'yscale','log','ylim',[5e-11,2e-7],'xlim',[0 4],'xtick',[],'ytick',[1e-10,1e-8],'FontSize',15)
end
for i=4:6
    subplot(2,3,i)
    set(gca,'yscale','log','ylim',[5e-11,1e-5],'xlim',[0 4],'xticklabel',{'0,18','0,300','0,18,300'},'FontSize',15)
end

figure(12); clf;
AA = [detCovarTotal(1,:),detCovarTotal(4,:),detCovarTotal(5,:)];
BB = [1./[Initial_Total],1./[Initial_MCP],1./[Initial_smiFISH]];
scatter(AA,BB,200,'filled')
set(gca,'YScale','log','xscale','log')

figure(13); clf;
mkr = ['o';'s';'d';'v'];
col = ['k','m','c'];
for i=1:4
    for j=1:3
        subplot(1,3,1)
        plot(InitialsigsTotal((j-1)*4+i),sigsTotal(1,(j-1)*4+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(j)); hold on
        set(gca,'YScale','log','xscale','log')
        subplot(1,3,2)
        plot(InitialsigsMCP((j-1)*4+i),sigsTotal(4,(j-1)*4+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(j)); hold on
        set(gca,'YScale','log','xscale','log')
        subplot(1,3,3)
        plot(InitialsigsFISH((j-1)*4+i),sigsTotal(5,(j-1)*4+i),[mkr(i)],'markersize',16,'MarkerFaceColor',col(j)); hold on
        set(gca,'YScale','log','xscale','log')
    end
end

%% Make Final MH Uncertainty plots
close all
vars.doFit = 0;
vars.priorScale = 1;
vars.modelVarsToFit = [2:5];
vars.mleScatterPlots = false;
FittingFunctionsCoLocalized(4,[FN,'_0_300'],vars);
FittingFunctionsCoLocalized(4,[FN,'_0_18_300'],vars);
bestFitTotal = FittingFunctionsCoLocalized(4,[FN,'_0_18_300'],vars);

for iFig = 2:3:14
    figure(iFig);
    for i =[1,2,3,5,6,9]
        h=subplot(3,3,i);
        axi = mod(i-1,3)+2;
        axj = floor((i-1)/4+1);
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');
        plot((bestFitTotal.MLE(1,axi))*[1,1],ylim+[-3,3])
        plot(xlim+[-3,3],(bestFitTotal.MLE(1,axj))*[1,1])
    end
end

%%         Reformat Scatter Plots
MLD_lines_to_plot = [2,7,3,1]+2;
MLD_lines_to_skip = [7,8];
for iFig = [2,11,14]
    figure(iFig);
    for i =[1,2,3,5,6,9]
        h=subplot(3,3,i);
        for j=MLD_lines_to_plot
            h.Children(j).Visible = 'on';
            h.Children(j).set('LineWidth',3)
        end
        for j=MLD_lines_to_skip
            h.Children(j).Visible = 'off';
        end
        h.Children(1).set("Color",[.6,.6,.6],'LineStyle','--')
        h.Children(2).set("Color",[.6,.6,.6],'LineStyle','--')
        h.Children(3).set("Color",'k',"Marker",'diamond','LineStyle','none')
        h.Children(4).set("Color",'c')
        h.Children(5).set("Color",'k','LineStyle','-')
        h.Children(6).set("Color",'k','LineStyle','--')
         h.Children(9).set("Color",'m','LineStyle','-')
         h.Children(10).set("Color",'m','LineStyle','--')
    
        if i==16
        legend(h.Children(MLD_lines_to_plot),...
            {'Prior','\Sigma, t=\{0\}',...
            '\Sigma, t=\{0,18\}','\Sigma, t=\{0,300\}',...
            '\Sigma, t=\{0,18,300\}','MLE, t=\{0,18,300\}'})
        end  
        set(gca,'FontSize',16)
    end

ax = [-4.1,-3.4; -2,-0.3; 1, 2.2; -2.8,-1.4];
for i=1:3
    subplot(3,3,i)
    set(gca,'ylim',ax(1,:))
end
subplot(3,3,1)
set(gca,'xlim',ax(2,:))
for i=5:6
    subplot(3,3,i)
    set(gca,'ylim',ax(2,:))
end
for i=3:3:9
    subplot(3,3,i)
    set(gca,'xlim',ax(4,:))
end
for i=2:3:5
    subplot(3,3,i)
    set(gca,'xlim',ax(3,:))
end
subplot(3,3,9)
set(gca,'ylim',ax(3,:))
end

 
%% Make Tex Table for Final parameters;
load([FN,'_0_18_300']);
Atxt = importdata('InitParTableTex.txt');
va = ['O';'W';'B';'G'];
vn = [2:5];
pdo = [1,4,5];

for i=1:4
    for j=1:3
        Atxt = strrep(Atxt,[va(i),'1',num2str(j)],num2str(log10(ModelZero{pdo(j)}.parameters{vn(i),2}),'%0.3f'));

        smplDone = chainResults{pdo(j)}.mhSamples(:,:);
        valDone = chainResults{pdo(j)}.mhValue;
        smplDone = smplDone(valDone~=0,:);
        covMH = cov(smplDone/log(10));
        Atxt = strrep(Atxt,[va(i),'2',num2str(j)],num2str(sqrt(covMH(i,i)),'%0.3f'));

        covFIM = covLog{pdo(j)}/(log(10)^2);
        Atxt = strrep(Atxt,[va(i),'3',num2str(j)],num2str(sqrt(covFIM(i,i)),'%0.3f'));

    end
end
Atxt

%% Functions

function [Final_Total,Sigs_Total] = evaluateFIMDesigns(FIMR_Final,iPDO,NCells,pars)
for i=1:size(FIMR_Final.fimResultsOne,1)
    FIM = NCells(1)*FIMR_Final.fimResultsOne{i,iPDO}{2};
    FIM = diag(FIMR_Final.parsets{iPDO}(i,:))*FIM*diag(FIMR_Final.parsets{iPDO}(i,:));
    detFIM_Total_Best_0(i,:) = det(FIM(pars,pars));
    sigsFIM_Total_Best_0(i,:) = diag(FIM(pars,pars)^(-1));
    

    FIM = NCells(1)*FIMR_Final.fimResultsOne{i,iPDO}{2}+...
        NCells(2)*FIMR_Final.fimResultsOne{i,iPDO}{6};
    FIM = diag(FIMR_Final.parsets{iPDO}(i,:))*FIM*diag(FIMR_Final.parsets{iPDO}(i,:));
    detFIM_Total_Best_0_18(i,:) = det(FIM(pars,pars));
    sigsFIM_Total_Best_0_18(i,:) = diag(FIM(pars,pars)^(-1));

    FIM = NCells(1)*FIMR_Final.fimResultsOne{i,iPDO}{2}+...
        NCells(3)*FIMR_Final.fimResultsOne{i,iPDO}{53};
    FIM = diag(FIMR_Final.parsets{iPDO}(i,:))*FIM*diag(FIMR_Final.parsets{iPDO}(i,:));
    detFIM_Total_Best_0_300(i,:) = det(FIM(pars,pars));
    sigsFIM_Total_Best_0_300(i,:) = diag(FIM(pars,pars)^(-1));

    FIM = NCells(1)*FIMR_Final.fimResultsOne{i,iPDO}{2}+...
        NCells(2)*FIMR_Final.fimResultsOne{i,iPDO}{6}+...
        NCells(3)*FIMR_Final.fimResultsOne{i,iPDO}{53};
    FIM = diag(FIMR_Final.parsets{iPDO}(i,:))*FIM*diag(FIMR_Final.parsets{iPDO}(i,:));
    detFIM_Total_Best_0_18_300(i,:) = det(FIM(pars,pars));
    sigsFIM_Total_Best_0_18_300(i,:) = diag(FIM(pars,pars)^(-1));
end

% Final_Total = [mean((detFIM_Total_Best_0)),mean((detFIM_Total_Best_0_18)),...
%     mean((detFIM_Total_Best_0_300)),mean((detFIM_Total_Best_0_18_300))];
Final_Total = [mean((detFIM_Total_Best_0_18)),...
    mean((detFIM_Total_Best_0_300)),mean((detFIM_Total_Best_0_18_300))];
Sigs_Total = [mean((sigsFIM_Total_Best_0_18)),...
    mean((sigsFIM_Total_Best_0_300)),mean((sigsFIM_Total_Best_0_18_300))];
end
