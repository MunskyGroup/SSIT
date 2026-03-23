%% example_SI_Epidemics
% Example script to demonstrate modeling epidemic data using SI, SIS, SIR, 
% and SEIS epidemiological models
%clear; clc; close all
addpath(genpath('../src'));

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

% Set solution scheme to 'ODE':
SI.solutionScheme = 'ODE';
SIS.solutionScheme = 'ODE';
SIR.solutionScheme = 'ODE';
SEIR.solutionScheme = 'ODE';
    
% Solve ODEs    
SI.Solutions = SI.solve; 
SIS.Solutions = SIS.solve; 
SIR.Solutions = SIR.solve; 
SEIR.Solutions = SEIR.solve; 

% Plot ODE solutions
SI.plotODE(speciesNames=SI.species, timeVec=SI.tSpan)
SIS.plotODE(speciesNames=SIS.species, timeVec=SIS.tSpan)
SIR.plotODE(speciesNames=SIR.species, timeVec=SIR.tSpan)
SEIR.plotODE(speciesNames=SEIR.species, timeVec=SEIR.tSpan)

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
[SI_FSPsoln,SI.fspOptions.bounds] = SI.solve; 
[SIS_FSPsoln,SIS.fspOptions.bounds] = SIS.solve; 
[SIR_FSPsoln,SIR.fspOptions.bounds] = SIR.solve;
[SEIR_FSPsoln,SEIR.fspOptions.bounds] = SEIR.solve;

% Plot marginal distributions
SI.plotFSP(solution=SI_FSPsoln, plotType='marginals', indTimes=[1:5:50])
SIS.plotFSP(solution=SIS_FSPsoln, plotType='marginals', indTimes=[1:5:50])
SIR.plotFSP(solution=SIR_FSPsoln, plotType='marginals', indTimes=[1:5:50])
SEIR.plotFSP(solution=SEIR_FSPsoln,plotType='marginals',indTimes=[1:5:50])
