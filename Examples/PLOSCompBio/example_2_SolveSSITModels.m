%% example_SolveSSITModels
% Example script to show how to solve the time evolution of state space 
% probabilities for a reaction system:  
% * Deterministically using ordinary differential equations (ODEs)
% * Stochastically using stochastic simulation algorithm (SSA) trajectories 
% * and Finite State Projection (FSP) of the Chemical Master Equation (CME)
% * With a hybrid approach, using upstream ODEs and downstream FSP.
clear
close all
addpath(genpath('../../'));

%% Preliminary
% Load our models from example_1_CreateSSITModels and inspect them
example_1_CreateSSITModels

Model.summarizeModel
STL1_Model.summarizeModel

%% Ex.(1) Compute Ordinary Differential Equations (ODEs)
% Create a copy of the Model
Model_ODE = Model;

% Set solution scheme to 'ode'
Model_ODE.solutionScheme = 'ode';


Model_ODE = Model_ODE.formPropensitiesGeneral('Model_ODE');

% Solve ODE and make plots
ODEsoln = Model_ODE.solve; 
plotODE(ODEsoln,Model_ODE.species)

%% Ex.(2): Run Gillepsie's Stochastic Simulation Algorithm (SSA) 
%% and analyse trajectories


%% This model consists of 3 species: 
% an inactive gene ('offGene'), an activated gene ('onGene'), and mRNA. 

%% There are four reactions, each with a unique rate parameter: 
% activation ('kon'), inactivation ('koff'), transcription ('kr'), 
% and mRNA degradation ('gr')


% Set initial condition (one 'offGene'):
Model.initialCondition = [1;0;0]; 

% Print a summary of our Model:
Model.summarizeModel


%% Ex.(3) Use the Finite State Projection (FSP) approximation of the CME
%% Model:

Model = Model.formPropensitiesGeneral('Model'); % Generate model codes

% Select FSP solution Scheme
Model.solutionScheme = 'FSP';  

% Set FSP 1-norm error tolerance.
Model.fspOptions.fspTol = 1e-4; 

% Guess initial bounds on FSP StateSpace
Model.fspOptions.bounds(4:6) = [2,2,400];

% Set times (s) at which to compute distributions:
Model.tSpan = linspace(0,180,301);

Model.fspOptions.initApproxSS = false; 

% Solve Model
[FSPsoln,Model.fspOptions.bounds] = Model.solve; 

% Plot marginal distributions
Model.makePlot(FSPsoln,'marginals',[1:100:301],false,[1,2,3],{'linewidth',2})  
Model.makePlot(FSPsoln,'margmovie',[],false,[101],{'linewidth',2},'movie.mp4',[1,1,0.015],[2,3])  

%% STL1_Model:
STL1_Model = STL1_Model.formPropensitiesGeneral('STL1_Model'); % Generate model codes

% Select FSP solution Scheme
STL1_Model.solutionScheme = 'FSP';  

% Set FSP 1-norm error tolerance.
STL1_Model.fspOptions.fspTol = 1e-4; 

% Guess initial bounds on FSP StateSpace
STL1_Model.fspOptions.bounds(4:6) = [2,2,400];

% Set times (s) at which to compute distributions:
STL1_Model.tSpan = (0:10:100);

STL1_Model.fspOptions.initApproxSS = false; 

% Solve Model
[STL1_FSPsoln,STL1_Model.fspOptions.bounds] = STL1_Model.solve; 

% Plot marginal distributions
STL1_Model.makePlot(STL1_FSPsoln,'marginals',[1:100:100],false,[1,2,3],{'linewidth',2})  
STL1_Model.makePlot(STL1_FSPsoln,'margmovie',[],false,[101],{'linewidth',2},'movie.mp4',[1,1,0.015],[2,3])  

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
% Try tinkering with the MAPK signal parameters (parameters 5-8)to get it 
% to match somewhat qualitatively to what we see in experiments: maximum at
% ~2 minutes and adaptation in ~10 min.
% We don't have to be exact, ballpark parameters should be fine to start.

%% Solve and plot using the FSP approach
% To solve the model, we first select the solution scheme ('FSP') and then
% we call the SSIT.solve method.
Model.parameters = ({'k21',30;'kr',100;'deg',0.005; ...
    'a0',0.01;'a1',1;'r1',0.4;'r2',.1});

Model.solutionScheme = 'FSP';    % Set solutions scheme to FSP.

% Set the code to start at steady state at t=0;
Model.fspOptions.initApproxSS =true;
Model = Model.formPropensitiesGeneral('STL1Model',true);
[FSPsoln,Model.fspOptions.bounds] = Model.solve;  % Solve the FSP analysis

% Next we make plots of the marginal distributions at time points 3, 5, 7,
% 9, 11, 13 and plot these in figures 1:3 for the three different species.
Model.makePlot(FSPsoln,'marginals',[3:2:13],false,(1:3))    % Plot marginal distributions

% We can also plot the means and standard deviations versus time in figure
% 100:
Model.makePlot(FSPsoln,'meansAndDevs',[],false,100)    % Plot marginal distributions

% Try to tune the parameters until you see:
% Bimodal expression (i.e., a population of active cells and a population of
% inactive cells).
% Perfect adaptation (all mRNA gone) at about 25 min.
% An average of ~50 mRNA at the highest expression time.