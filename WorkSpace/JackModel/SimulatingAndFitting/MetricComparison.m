%% Simulate diffusion
clear
close all 
addpath(genpath('../../../src'));
addpath("RNA_Modeling")


%% 
kon = 1;
koff = 1;
w = 1;
kex = 100;
kr = 100;
D = [0.01,10,10];
gam =[1;1;1];
Ncells = 100;
a = 53*ones(Ncells,1);
b = 39*ones(Ncells,1);
u = 2*pi*rand(Ncells,1);
posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data
makePlot = 1;
fileName = 'regime1.csv';

stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

[intensTS_R1,binCounts_R1,~,~,~,binTS_R1,binCounts_R1Dim,binCounts_R1Bright,...
    hisCountsBright_R1,hisCountsDim_R1,XcRot_R1,iTS_R1,XCVcRot_R1,XYfit_R1] = ...
    processSpotPositionData(fileName,[],a,b,0);

kon = 1;
koff = 1;
w = 1;
kex = 100;
kr = 100;
D = [0.01,10,10];
gam =[1;1;1];
Ncells = 100;
a = 53*ones(Ncells,1);
b = 39*ones(Ncells,1);
u = 2*pi*rand(Ncells,1);
posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data
makePlot = 1;
fileName = 'regime2.csv';

stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

[intensTS_R2,binCounts_R2,~,~,~,binTS_R2,binCounts_R2Dim,binCounts_R2Bright,...
    hisCountsBright_R2,hisCountsDim_R2,XcRot_R2,iTS_R2,XCVcRot_R2,XYfit_R2] = ...
    processSpotPositionData(fileName,[],a,b,0);

kon = 1;
koff = 1;
w = 1;
kex = 100;
kr = 100;
D = [0.01,10,10];
gam =[1;1;1];
Ncells = 100;
a = 53*ones(Ncells,1);
b = 39*ones(Ncells,1);
u = 2*pi*rand(Ncells,1);
posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data
makePlot = 1;
fileName = 'regime1.csv';

stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

[intensTS_R3,binCounts_R3,~,~,~,binTS_R3,binCounts_R3Dim,binCounts_R3Bright,...
    hisCountsBright_R3,hisCountsDim_R3,XcRot_R3,iTS_R3,XCVcRot_R3,XYfit_R3] = ...
    processSpotPositionData(fileName,[],a,b,0);

kon = 1;
koff = 1;
w = 1;
kex = 100;
kr = 100;
D = [0.01,10,10];
gam =[1;1;1];
Ncells = 100;
a = 53*ones(Ncells,1);
b = 39*ones(Ncells,1);
u = 2*pi*rand(Ncells,1);
posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data
makePlot = 1;
fileName = 'regime1.csv';

stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

[intensTS_R4,binCounts_R4,~,~,~,binTS_R4,binCounts_R4Dim,binCounts_R4Bright,...
    hisCountsBright_R4,hisCountsDim_R4,XcRot_R4,iTS_R4,XCVcRot_R4,XYfit_R4] = ...
    processSpotPositionData(fileName,[],a,b,0);



kon = 1;
koff = 1;
w = 1;
kex = 100;
kr = 100;
D = [0.01,10,10];
gam =[1;1;1];
Ncells = 100;
a = 53*ones(Ncells,1);
b = 39*ones(Ncells,1);
u = 2*pi*rand(Ncells,1);
posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data
makePlot = 1;
fileName = 'regime1.csv';

stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

[intensTS_R5,binCounts_R5,~,~,~,binTS_R5,binCounts_R5Dim,binCounts_R5Bright,...
    hisCountsBright_R5,hisCountsDim_R5,XcRot_R5,iTS_R5,XCVcRot_R5,XYfit_R5] = ...
    processSpotPositionData(fileName,[],a,b,0);