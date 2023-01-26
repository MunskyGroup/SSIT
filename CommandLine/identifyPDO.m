close all

DATA = importdata('LindaDataAllTimes.xlsx');
I = find(strcmp(DATA.colheaders,'RNA_nuc_ch1_ts500'));
J = find(strcmp(DATA.colheaders,'RNA_nuc_ch3_ts500'));

YI = DATA.data(DATA.data(:,7)==0,J);
XI = DATA.data(DATA.data(:,7)==0,I);

OBJ = @(l)-findError(l,XI,YI);
options = optimset('display','iter');
LL12 = [.7;10;0.1];
for i=1:1
LL12 = fminsearch(OBJ,LL12,options);
end

OBJ = @(l)-findError(l,YI,XI);
options = optimset('display','iter');
LL21 = [.7;10;0.1];
for i=1:1
LL21 = fminsearch(OBJ,LL21,options);
end

%%
close all
[logL,P12] = findError(LL12,XI,YI);
[logL,P21] = findError(LL21,YI,XI);

figure(1)
contourf(log(P12),-450:50:-50)
colorbar
hold on
scatter(XI,YI,100,'sm','filled')
xlabel('Number observed in Cond 1')
ylabel('Number observed in Cond 2')
plot([0,min(size(P12))],[0,min(size(P12))],'k--','LineWidth',3)
set(gca,'fontsize',16)

figure(2)
contourf(log(P21))
colorbar
hold on
scatter(YI,XI,100,'sm','filled')
xlabel('Number observed in Cond 1')
ylabel('Number observed in Cond 2')
plot([0,min(size(P12))],[0,min(size(P12))],'k--','LineWidth',3)
set(gca,'fontsize',16)

%%
col1 = [0.2,0.2,0.6];
col2 = [0.6,0.2,0.2];
figure(3); clf; 
Nmax = min(400,size(P12,2));
edges = 0:Nmax;
histogram(XI,edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col1); hold on
histogram(YI,edges,'normalization','cdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col2);
set(gca,'FontSize',16,'xlim',[0,Nmax-1])

figure(4); clf
binSize = 20;
Nmax = min(400,size(P12,2));
edges = 0:binSize:Nmax;
histogram(XI,edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col1); hold on
histogram(YI,edges,'normalization','pdf','DisplayStyle','stairs','LineWidth',4,'EdgeColor',col2);
set(gca,'FontSize',16,'xlim',[0,Nmax-1])

figure(10); clf
edges = 0:Nmax;
h1 = histogram(XI,edges,'normalization','pdf');
h1 = h1.Values;
h2 = histogram(YI,edges,'normalization','pdf');
h2 = h2.Values;
close(10);

C12 = P12(1:Nmax,1:Nmax);
C21 = P21(1:Nmax,1:Nmax);
h2Pred = C12*h1';
h1Pred = C21*h2';

P2modBin=[];
P1modBin=[];
for i = 1:floor(Nmax/binSize)
    P2modBin(i) = sum(h2Pred((i-1)*binSize+1:i*binSize))/binSize;
    P1modBin(i) = sum(h1Pred((i-1)*binSize+1:i*binSize))/binSize;
end
figure(3)
stairs([0:Nmax],[cumsum(h2Pred);1],'--','LineWidth',4,'Color',col2);
stairs([0:Nmax],[cumsum(h1Pred);1],'--','LineWidth',4,'Color',col1);


figure(4)
stairs([0:binSize:Nmax-binSize],(P2modBin),'--','LineWidth',4,'Color',col2);
stairs([0:binSize:Nmax-binSize],(P1modBin),'--','LineWidth',4,'Color',col1);



% % figure(2); clf
% S1 = [smooth(H1.Values,5)',0];
% S2 = [smooth(H2.Values,5)',0];
% stairs(edges,S1,'m','LineWidth',4); hold on
% stairs(edges,S2,'color',[0.2 0.5 0.2],'LineWidth',4); hold on
% legend('smFISH','MCP-GFP')
% set(gca,'FontSize',16)
% xlabel('Number of mRNA')
% ylabel('Empirical PMF')

%%
function [logL,P] = findError(lambda,XI,YI)
% Computes likelihood of observed data given the model of binomial spot
% loss and poisson extra spot counting.
Nmax = max(max(XI,YI));
Np = max(max(YI-XI)+50,lambda(2)+ceil(5*sqrt(lambda(2))));
P = zeros(Nmax+Np,Nmax);
for xi = 0:Nmax
    P2 = pdf('poiss',[0:Np],lambda(2)+lambda(3)*xi);
    P1 = binopdf([0:xi],xi,lambda(1));
    P(1:xi+Np+1,xi+1) = conv(P2,P1);
end

%%
logP = log(P);

logL = 0;
for i = 1:length(XI)
    logL = logL + logP(YI(i)+1,XI(i)+1);
end
logL = logL-1e4*((lambda(1)>1)+(lambda(1)<0)+(lambda(2)<0));
end

% %%
% mnY = mean(YI)
% w = 0.2;
% b = 15.4*0.67/(0.67+0.78);
% k21range = -log(0.5)./[100,400];
% trange = [30,180];
% grange = -log(0.5)./trange;
% k12range = mn
% 
% varY = var(YI)
% 

