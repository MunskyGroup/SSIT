Forman = Forman2027;
%% Figure 1 B, C, D 
f1B=figure(101); ax1B = gca;
Forman.makeFig1B(ax1B);

f1C1=figure(102);
f1C2=figure(103);
f1C3=figure(104); 
f1D=figure(105);
Forman.makeFigs1CD(f1C1,f1C2,f1C3,f1D);

% Forman.exportFigs1(f1B,f1C1,f1C2,f1C3,f1D)

%% Figure 2 A-H
figure(201); ax2A = gca;
Forman.makeFigs2A(ax2A);

Forman.prepareFigs2B;
%%
figure(202); ax2B = gca;
figure(203); ax2C = gca;
figure(204); ax2D = gca;
f2E1 = figure(205);
figure(206); ax2E2 = gca;
figure(207); ax2F1 = gca;
figure(208); ax2F2 = gca;
Forman.makeFigs2BtoF(ax2B,ax2C,ax2D,f2E1,ax2E2,ax2F1,ax2F2,Nsamps=20);

%% 2G - Prepare MLE reuslts
Forman.prepareFig2G(nMLE=20)
%% 2G - Make the figures
Forman.makeFig2G;
Forman.makeFig2H;

%% Figure 3 A
f3A1 = figure(301); 
f3A2 = figure(302);
f3A3 = figure(303);
f3A4 = figure(304);
f3A5 = figure(305);

Forman.makeFigs3A(f3A1,f3A2,f3A3,f3A4,f3A5);

%% Compute FIMs for all possible experiment time points
Forman = Forman.prepareFIMS;

%% Make Figures 3B-D
f3B1 = figure(306); 
f3B2 = figure(307);
f3B3 = figure(308);
f3B4 = figure(309);
Forman = Forman.optimizeExperiment;
Forman = Forman.makeFig3B(f3B1,f3B2,f3B3,f3B4);
Forman = Forman.computeFIMDeterminants;

f3C1 = figure(310); 
f3C2 = figure(311);
f3C3 = figure(312);
f3C4 = figure(313);
Forman = Forman.makeFig3C(f3C1,f3C2,f3C3,f3C4);

f3D1 = figure(314); 
f3D2 = figure(315);
f3D3 = figure(316);
Forman = Forman.makeFig3D(f3D1,f3D2,f3D3);
fignums = [f3A1,f3A2,f3A3,f3A4,f3A5,f3B1,f3B2,f3B3,f3B4,f3C1,f3C2,f3C3,f3C4,f3D1,f3D2,f3D3];
Forman.exportFigs3(fignums,'AnnualReview_Figures',fileType='pdf')

%% Figure 4ABC 
f4A = figure(401);
f4B = figure(402);
f4C = figure(403);
Forman = Forman.makeFig4ABC(f4A,f4B,f4C);

%% Figure 4DEF
Forman = Forman.prepareFig4DEF(nMLE=200);
%%
Forman.makeFigs4DEF;