Forman = Forman2027;
%% Figure 1 B, C, D 
f1B=figure(101); ax1B = gca;
Forman.makeFig1B(ax1B);

%%
f1C1=figure(102);
f1C2=figure(103);
f1C3=figure(104); 
f1D=figure(105);
Forman.makeFigs1CD(f1C1,f1C2,f1C3,f1D);

% Forman.exportFigs1(f1B,f1C1,f1C2,f1C3,f1D)

%% Figure 2 A-H
f2a = figure(201); ax2A = gca;
Forman.makeFigs2A(ax2A);

Forman.prepareFigs2B;
%%
f2b = figure(202); ax2B = gca;
f2c = figure(203); ax2C = gca;
f2d = figure(204); ax2D = gca;
f2E1 = figure(205);
f2e = figure(206); ax2E2 = gca;
figure(207); ax2F1 = gca;
f2f = figure(208); ax2F2 = gca;
Forman.makeFigs2BtoF(ax2B,ax2C,ax2D,f2E1,ax2E2,ax2F1,ax2F2,Nsamps=20);

%% 2G - Prepare MLE reuslts
Forman.prepareFig2G(nMLE=500)
%% 2G - Make the figures
f2g1 = figure(209);
f2g2 = figure(210);
f2h1 = figure(211);
f2h2 = figure(212);
Forman.makeFig2G;
Forman.makeFig2H(f2g1, f2g2, f2h1, f2h2);

%% 
fignums = [f2a, f2b, f2c, f2d, f2e, f2f, f2g1, f2g2, f2h1, f2h2];
Forman.exportFigs3(fignums,'AnnualReview_Figures',fileType='svg')
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
Forman.exportFigs3(fignums,'AnnualReview_Figures',fileType='svg')

%% Figure 4ABC 
f4A = figure(401);
f4B = figure(402);
f4C = figure(403);
Forman = Forman.makeFig4ABC(f4A,f4B,f4C);


%% Figure 4DEF
Forman.freeParsFig4 = [1:4];
Forman = Forman.prepareFig4DEF(nMLE=250);


%%
f4d1 = figure(404);
f4d2 = figure(405);
f4e1 = figure(406);
f4e2 = figure(407);
f4f1 = figure(408);
f4f2 = figure(409);
Forman.makeFigs4DEF(f4d1, f4d2, f4e1, f4e2, f4f1, f4f2);


%% 
f4g = figure(410);
f4h = figure(411);
f4i = figure(412);
Forman.makeFig4GHI(f4g, f4h, f4i);


%%
% fignums = [f4A,f4B,f4C,f4d1,f4d2,f4e1,f4e2,f4f1,f4f2, f4g, f4h, f4i];
fignums = [f4g, f4h, f4i];
Forman.exportFigs3(fignums,'AnnualReview_Figures',fileType='svg')
