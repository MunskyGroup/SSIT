%% Regime Exploration
clear
close all 
addpath(genpath('../../../src'));

%% Pretext 
% Iyer-Biswas, Srividya, et al. “Stochasticity of Gene Products from Transcriptional Pulsing.”
% This paper explored the parameter space of the one compartment model,
% some key finding was that you can fix everything relative to one
% constant (kd) and explore the parameters, kr/kd controls the extent of the
% distribution and not the shape, kon and koff primarly control the shape
% we will explore the five regimes:
% regime 1 - kd = kon = koff
% regime 2 - kon, koff > kd
% regime 3 - kon, koff < kd
% regime 4 - kon < kd < koff
% regime 5 - koff < kd < kon

% this file will explore how two compartements effects the shape of these
% regimes as a function of the transition rate

%% Build Model 
% One Compartment Model
OneCompartmentModel = SSIT();
OneCompartmentModel.parameters = {'kon',1; 'koff',1; 'kr',100; 'kd',1};
OneCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'};
OneCompartmentModel.stoichiometry = [1, -1, 0, 0;...
                                    -1, 1, 0, 0;...
                                    0, 0, 1, -1];
OneCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'kd*RNA'};
OneCompartmentModel.initialCondition = [0;1;0;];
OneCompartmentModel.summarizeModel
OneCompartmentModel = OneCompartmentModel.formPropensitiesGeneral('OneComparment_SpatialModel');
OneCompartmentModel.tSpan = linspace(0,10,21);
% OneCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on';'gene_off';};
save('OneComparment_SpatialModel','OneCompartmentModel')

% Two Compartment Model
TwoCompartmentModel = SSIT();
TwoCompartmentModel.parameters = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
TwoCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'};
TwoCompartmentModel.stoichiometry = [1, -1, 0, 0, 0, 0, 0;... % species x reaction
                                    -1, 1, 0, 0, 0, 0, 0;...
                                    0, 0, 1, -1, 0, -1, 1;...
                                    0, 0, 0, 1, -1, 0, -1;];
TwoCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'D*RNA'; 'kd*D_RNA'; 'kd*RNA'; 'D*D_RNA';};
TwoCompartmentModel.initialCondition = [0;1;0;0;];
TwoCompartmentModel.summarizeModel
TwoCompartmentModel = TwoCompartmentModel.formPropensitiesGeneral('TwoComparment_SpatialModel');
TwoCompartmentModel.tSpan = linspace(0,10,21);
% TwoCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on';'gene_off';};
save('TwoComparment_SpatialModel','TwoCompartmentModel')



%% Regime 1
% 1 compartment
model1 = SSIT('OneComparment_SpatialModel','OneCompartmentModel');
model1.parameters = {'kon',1; 'koff',1; 'kr',100; 'kd',1};
model1.solutionScheme = 'FSP';
[FSPsoln,model1.fspOptions.bounds] = model1.solve; 
fig1 = figure(101);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(102);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(103);clf; set(fig3,'Name','Marginal Distributions, RNA');
model1.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3],{'r-','linewidth',2})

% D = 1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(111);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(112);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(113);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(114);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})

% D = 10
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1;'koff',1;'kr',100; 'D', 10; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(121);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(122);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(123);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(124);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})


% D = 0.1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1;'koff',1;'kr',100; 'D', 0.1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(131);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(132);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(133);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(134);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})



%% Regime 2
% 1 compartment
model1 = SSIT('OneComparment_SpatialModel','OneCompartmentModel');
model1.parameters = {'kon',1.2; 'koff',1.5; 'kr',100; 'kd',1};
model1.solutionScheme = 'FSP';
[FSPsoln,model1.fspOptions.bounds] = model1.solve; 
fig1 = figure(201);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(202);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(203);clf; set(fig3,'Name','Marginal Distributions, RNA');
model1.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3],{'r-','linewidth',2})

% D = 1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1.2;'koff', 1.5;'kr',100; 'D', 1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(211);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(212);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(213);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(214);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})

% D = 10
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1.2;'koff', 1.5;'kr',100; 'D', 10; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(221);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(222);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(223);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(224);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})


% D = 0.1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1.2;'koff', 1.5;'kr',100; 'D', 0.1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(231);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(232);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(233);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(234);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})


%% Regime 3
% 1 compartment
model1 = SSIT('OneComparment_SpatialModel','OneCompartmentModel');
model1.parameters = {'kon',0.5; 'koff',0.5; 'kr',100; 'kd',1};
model1.solutionScheme = 'FSP';
[FSPsoln,model1.fspOptions.bounds] = model1.solve; 
fig1 = figure(301);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(302);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(303);clf; set(fig3,'Name','Marginal Distributions, RNA');
model1.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3],{'r-','linewidth',2})

% D = 1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',0.5;'koff', 0.5;'kr',100; 'D', 1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(311);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(312);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(313);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(314);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})

% D = 10
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',0.5;'koff', 0.5;'kr',100; 'D', 10; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(321);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(322);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(323);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(324);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})


% D = 0.1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',0.5;'koff', 0.5;'kr',100; 'D', 0.1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(331);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(332);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(333);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(334);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})



%% Regime 4
% 1 compartment
model1 = SSIT('OneComparment_SpatialModel','OneCompartmentModel');
model1.parameters = {'kon',0.5; 'koff',2; 'kr',100; 'kd',1};
model1.solutionScheme = 'FSP';
[FSPsoln,model1.fspOptions.bounds] = model1.solve; 
fig1 = figure(401);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(402);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(403);clf; set(fig3,'Name','Marginal Distributions, RNA');
model1.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3],{'r-','linewidth',2})

% D = 1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',0.5;'koff', 2;'kr',100; 'D', 1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(411);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(412);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(413);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(414);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})

% D = 10
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',0.5;'koff', 2;'kr',100; 'D', 10; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(421);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(422);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(423);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(424);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})


% D = 0.1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',0.5;'koff', 2;'kr',100; 'D', 0.1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(431);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(432);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(433);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(434);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})


%% Regime 5
% 1 compartment
model1 = SSIT('OneComparment_SpatialModel','OneCompartmentModel');
model1.parameters = {'kon',1.5; 'koff',0.5; 'kr',100; 'kd',1};
model1.solutionScheme = 'FSP';
[FSPsoln,model1.fspOptions.bounds] = model1.solve; 
fig1 = figure(501);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(502);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(503);clf; set(fig3,'Name','Marginal Distributions, RNA');
model1.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3],{'r-','linewidth',2})

% D = 1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1.5;'koff', 0.5;'kr',100; 'D', 1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(511);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(512);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(513);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(514);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})

% D = 10
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1.5;'koff', 0.5;'kr',100; 'D', 10; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(521);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(522);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(523);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(524);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})


% D = 0.1
model2 = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel');
model2.parameters = {'kon',1.5;'koff', 0.5;'kr',100; 'D', 0.1; 'kd', 1};
model2.solutionScheme = 'FSP';
[FSPsoln,model2.fspOptions.bounds] = model2.solve; 
fig1 = figure(531);clf; set(fig1,'Name','Marginal Distributions, gene_on');
fig2 = figure(532);clf; set(fig2,'Name','Marginal Distributions, gene_off');
fig3 = figure(533);clf; set(fig3,'Name','Marginal Distributions, RNA');
fig4 = figure(534);clf; set(fig3,'Name','Marginal Distributions, D_RNA');
model2.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig1, fig2, fig3, fig4],{'r-','linewidth',2})



%% Conlusions 
% I probably need to look at more intermediate values for D but it seems
% that the shape can be effect by D but more so the extent of both RNA and
% D_RNA.







