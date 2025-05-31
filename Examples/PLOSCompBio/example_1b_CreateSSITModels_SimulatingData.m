%% example_1b_CreateSSITModels_SimulatingData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate data for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries:
%clear
%close all
addpath(genpath('../../'));

%example_1_CreateSSITModels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create a complicated model to simulate noisy real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a copy of the STL1 model for simulation:
STL1_sim = STL1;

% Update propensity function for the gene activation reaction:
STL1_sim.propensityFunctions{1} = 'offGene * IHog';

% Define the time-varying TF/MAPK input signal:
STL1_sim.inputExpressions = {'IHog',...
               'a0 + a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0) + sigma'};

% Add the new parameters from the TF/MAPK input signal:
STL1_sim.parameters = ({'koff',0.01; 'kr',100; 'gr',0.5; 'kleak',0.5; ...
            'a0',0.01; 'a1',10; 'r1',0.004; 'r2',.01; 'sigma',randn(1)*5});

% Set stoichiometry of reactions:
STL1_sim.stoichiometry = [-1,1,0,0,0;...
                           1,-1,0,0,0;...
                           0,0,1,-1,1]; 

% Set propensity functions:
STL1_sim.propensityFunctions{5} = 'kleak'; 

% Print a summary of STL1 Model:
STL1_sim.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate simulated time-varying STL1 yeast data for fitting later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of independent data sets to generate:
STL1_sim.ssaOptions.Nexp = 2;  

% Number of cells to include at each time point for each data set:
STL1_sim.ssaOptions.nSimsPerExpt = 200;

% Ensure the solution scheme is set to FSP (default):
    STL1_sim.solutionScheme = 'FSP';  

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1_FSP'
    STL1_sim = STL1_sim.formPropensitiesGeneral('STL1_sim_FSP');
    
    % Set FSP 1-norm error tolerance:
    STL1_sim.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1_sim.fspOptions.bounds = [2,2,400];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1_sim.fspOptions.initApproxSS = false; 
    
    % Solve Model:
    [STL1_sim_FSPsoln,STL1_sim.fspOptions.bounds] = STL1_sim.solve;

% Generate and save data:
dataTable = STL1_sim.sampleDataFromFSP(STL1_sim_FSPsoln,'data/STL1_sim.csv'); 

% Plot data as histograms:
for i = 1:4
    subplot(2,2,i)  % Switch to current subplot
    histogram(dataTable.exp1_s3(dataTable.time==STL1.tSpan(i)),30,...
                                  "DisplayStyle","stairs")
end 