clear
close all
addpath(genpath('../../src'));

%% (1) Create a Model in SSIT  
Model = SSIT;  % Create an SSIT instance and give it a name. Here, 'Model'.

Model.species = {'offGene';'onGene';'mRNA'}; % Set species names.

% Set stoichiometry of reactions:
Model.stoichiometry = [-1,1,0,0;...
                        1,-1,0,0;...
                        0,0,1,-1]; 

% Define a time-varying TF/MAPK input signal:
Model.inputExpressions = {'IHog','(a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0))'};

% Set propensity functions:
Model.propensityFunctions = {'offGene*IHog';'k21*onGene';'kr*onGene';'deg*mRNA'}; 

% Set initial condition (one offgene):
Model.initialCondition = [1;0;0]; 

% Set times (s) at which to compute distributions:
Model.tSpan = [0:5:60]; 

%% Plot the TF/MAPK signal
% Next, we have to gues some initial guesses for parameters.
% First, let's tinker with the MAPK signal to get it to match somewhat
% qualitatively to what we see in experiments.  We don't have to be exact,
% ballpark parameters should be fine to start.
Model.parameters = ({'k21',30;'kr',100;'deg',0.005; ...
    'a0',0.01;'a1',1;'r1',0.4;'r2',.1});
par = [Model.parameters{:,2}];
t = [0:60];
TF = par(4)+par(5)*exp(-par(6)*t).*(1-exp(-par(7)*t)).*(t>0);
figure(1); plot(t,TF,'linewidth',3); 
set(gca,'fontsize',16)
xlabel('Time (min)'); ylabel('Hog1(t)')

%%      (1B) Extending a model with additional reactions
F3 = F1.addSpecies('protein');
F3 = F3.addParameter({'kpact',4;'gpact',1});
F3 = F3.addReaction('kpact*x1',[0;1]); 
F3 = F3.addReaction('gpact*protein',[0;-1]); 

%% (2) Solving a Model
%%      (2A) Solve using FSP
F1 = F1.formPropensitiesGeneral('BasicModel');
F1.tSpan = [-1:1:10];
F1.initialTime = -1;
F1.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
[FSPsoln,F1.fspOptions.bounds] = F1.solve;  % Solve the FSP analysis

%%          (2A.1) Solve FSP Model again using the bounds from the last solution
% If we start with the bounds computed in the first analysis, the solution is
% often much faster.
[FSPsoln] = F1.solve(FSPsoln.stateSpace);  % Solve the FSP analysis

%%          (2A.2) Make plots of FSP solution
F1.makePlot(FSPsoln,'meansAndDevs',[],[],1,{'linewidth',3,'color',[0,1,1]}) % Make plot of mean vs. time.
F1.makePlot(FSPsoln,'marginals',[],[],2,{'linewidth',3,'color',[0,0,1]}) % Make plot of mean vs. time.

%%      (2B) Solve using SSA
F2 = F1;
F2.solutionScheme = 'SSA';
SSASoln = F2.solve;

%%          (2B.1) Make plots of SSA solution
F2.makePlot(SSASoln,'trajectories',[],[],4) % Make some plots.
F1.makePlot(FSPsoln,'meansAndDevs',[],[],4,...
    {'linewidth',4,'color',[0,1,1],'Marker','s','MarkerSize',20}) % Add FSP Solution to plot.

%% (3) Sensitivity Analysis using FSP
F4 = F1;
F4.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = F4.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem

%%      (3A) Make plot of sensitivities
F4.makePlot(sensSoln,'marginals',[],[],4,{'linewidth',3,'color',[0,0,1]}) % Plot marginal sensitivities

