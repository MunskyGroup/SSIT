%% example_SI_Epidemics
% Example script to demonstrate modeling epidemic data using SI, SIS, SIR, 
% and SEIS epidemiological models
clear; clc; close all
addpath(genpath('../../src'));

%% Define SI Model
% Set up a simple model where susceptible (S) individuals become infected 
% (I) at some rate kI

% Create SI model
SI = SSIT;  
% Set species (effective population type) names
SI.species = {'S';'I'}; 
% Set initial condition
SI.initialCondition = [200;1];           
% Define propensity functions
SI.propensityFunctions = {'kI * S'}; 
% Define stoichiometry
SI.stoichiometry = [-1;1]; 
% Specify parameter guesses
SI.parameters = ({'kI',1});  
% Set Initial Distribution to Steady State
SI.fspOptions.initApproxSS = false;  
% Print visual summary of model
SI.summarizeModel     
% Save propensity functions
SI = SI.formPropensitiesGeneral('SI'); 
% Set times at which distributions will be computed
SI.tSpan = linspace(0,20,200);

%% Define SIS Model
% Set up a simple model where susceptible (S) individuals become infected 
% (I) at rate kI and become susceptible again at rate kS

% Create SIS model
SIS = SSIT;  
% Set species (effective population type) names
SIS.species = {'S';'I'}; 
% Set initial condition
SIS.initialCondition = [200;1];           
% Define propensity functions
SIS.propensityFunctions = {'kI*S';'kS*I'}; 
% Define stoichiometry
SIS.stoichiometry = [-1,1;1,-1]; 
% Specify parameter guesses
SIS.parameters = ({'kI',1;'kS',0.2});  
% Set Initial Distribution to Steady State
SIS.fspOptions.initApproxSS = false;  
% Print visual summary of model
SIS.summarizeModel     
% Save propensity functions
SIS = SIS.formPropensitiesGeneral('SIS'); 
% Set times at which distributions will be computed
SIS.tSpan = linspace(0,20,200);

%% Define SIR Model
% Set up a simple model where susceptible (S) individuals become infected 
% (I) at rate kI and are removed (R) - e.g., recovered with immunity, 
% quarantined, etc. - at rate kR

% Create SIR model
SIR = SSIT;  
% Set species (effective population type) names
SIR.species = {'S';'I';'R'}; 
% Set initial condition
SIR.initialCondition = [200;1;0];           
% Define propensity functions
SIR.propensityFunctions = {'kI*S';'kR*I'}; 
% Define stoichiometry
SIR.stoichiometry = [-1,0;
                     1,-1;
                     0,1]; 
% Specify parameter guesses
SIR.parameters = ({'kI',1;'kR',0.2});  
% Set Initial Distribution to Steady State
SIR.fspOptions.initApproxSS = false;  
% Print visual summary of model
SIR.summarizeModel     
% Save propensity functions
SIR = SIR.formPropensitiesGeneral('SIR'); 
% Set times at which distributions will be computed
SIR.tSpan = linspace(0,20,200);

%% Define SEIR Model
% Set up a simple model where susceptible (S) individuals become infected
% but are not yet infectious - i.e., enter a latency period (E) - at rate
% kE, become infectious (I) at rate kI, and are removed (R) - e.g., 
% recovered with immunity, quarantined, etc. - at rate kR

% Create SIR model
SEIR = SSIT;  
% Set species (effective population type) names
SEIR.species = {'S';'E';'I';'R'}; 
% Set initial condition
SEIR.initialCondition = [200;0;1;0];           
% Define propensity functions
SEIR.propensityFunctions = {'kE*S';'kI*E';'kR*I'}; 
% Define stoichiometry
SEIR.stoichiometry = [-1,0,0;
                     1,-1,0;
                     0,1,-1;
                     0,0,1]; 
% Specify parameter guesses
SEIR.parameters = ({'kE',0.5;'kI',0.5;'kR',0.2});  
% Set Initial Distribution to Steady State
SEIR.fspOptions.initApproxSS = false;  
% Print visual summary of model
SEIR.summarizeModel     
% Save propensity functions
SEIR = SEIR.formPropensitiesGeneral('SEIR'); 
% Set times at which distributions will be computed
SEIR.tSpan = linspace(0,20,200);

%% Compute Ordinary Differential Equations (ODEs)

% Set solution scheme to 'ode':
SI.solutionScheme = 'ode';
SIS.solutionScheme = 'ode';
SIR.solutionScheme = 'ode';
SEIR.solutionScheme = 'ode';
    
% Solve ODEs
SI_ODEsoln = SI.solve; 
SIS_ODEsoln = SIS.solve; 
SIR_ODEsoln = SIR.solve; 
SEIR_ODEsoln = SEIR.solve; 

% Plot ODE solutions
plotODE(SI_ODEsoln,SI.species,SI.tSpan)
plotODE(SIS_ODEsoln,SIS.species,SIS.tSpan)
plotODE(SIR_ODEsoln,SIR.species,SIR.tSpan)
plotODE(SEIR_ODEsoln,SEIR.species,SEIR.tSpan)

%% Solve CME using FSP
% Next, we can solve the model using the FSP.  In this example, we show how
% to run the code twice.  First call finds the FSP projection needed to
% solve the problem, and the second call solves using that projection.
% Select FSP solution Scheme
SI.solutionScheme = 'FSP';  
SIS.solutionScheme = 'FSP'; 
SIR.solutionScheme = 'FSP'; 
SEIR.solutionScheme = 'FSP'; 
% Set FSP 1-norm error tolerance
SI.fspOptions.fspTol = 1e-5; 
SIS.fspOptions.fspTol = 1e-5; 
SIR.fspOptions.fspTol = 1e-5; 
SEIR.fspOptions.fspTol = 1e-5; 
% Guess initial bounds on FSP StateSpace
SI.fspOptions.bounds(1:2) = [201,201]; 
SIS.fspOptions.bounds(1:3) = [201,201,201]; 
SIR.fspOptions.bounds(1:3) = [201,201,201]; 
SEIR.fspOptions.bounds(1:4) = [201,201,201,201];  
% Solve Model
[SIS_FSPsoln,SIS.fspOptions.bounds] = SIS.solve; 
[SIR_FSPsoln,SIR.fspOptions.bounds] = SIR.solve;
[SEIR_FSPsoln,SEIR.fspOptions.bounds] = SEIR.solve;

% Plot marginal distributions
SIS.makePlot(SIS_FSPsoln,'marginals',[1:20:200],false,[1,2],...
             {'linewidth',2})  
SIS.makePlot(SIS_FSPsoln,'margmovie',[],false,[101],{'linewidth',2},...
            'SIS.mp4',[1,1,0.015],[1,2]) 
SIR.makePlot(SIR_FSPsoln,'marginals',[1:20:200],false,[1,2,3],...
             {'linewidth',2})  
SIR.makePlot(SIR_FSPsoln,'margmovie',[],false,[101],{'linewidth',2},...
            'SIR.mp4',[1,1,0.015],[1,2,3]) 
SEIR.makePlot(SEIR_FSPsoln,'marginals',[1:20:200],false,[1,2,3,4],...
             {'linewidth',2})  
SEIR.makePlot(SEIR_FSPsoln,'margmovie',[],false,[101],{'linewidth',2},...
            'SEIR.mp4',[1,1,0.015],[1,2,3,4]) 

