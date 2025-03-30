%% example_5_FIMCalculation
% Example script to set up and solve the FSP-FIM matrix  
% with partial observations and probabilistic distortion.
clear
close all
addpath(genpath('../src'));

%% Preliminaries
% Load our models described in example_1_CreateSSITModels and  
% compute FSP solutions using example_2_SolveSSITModels_FSP
example_2_SolveSSITModels_FSP
Model.summarizeModel
STL1Model.summarizeModel

% Make copies of Model and STL1 Model
Model_sens = Model;
STL1_sens = STL1Model;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity
Model_sens.solutionScheme = 'fspSens'; 
STL1_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem
[sensSoln,bounds] = Model_sens.solve(FSPsoln.stateSpace); 
[STL1_sensSoln,STL1_bounds] = STL1_sens.solve(STL1_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis
% Model:
fig1 = figure(1);clf; set(fig1,'Name','Marginal Sensitivity, offGene');
fig2 = figure(2);clf; set(fig2,'Name','Marginal Sensitivity, onGene');
fig3 = figure(3);clf; set(fig3,'Name','Marginal Sensitivity, mRNA');
Model_sens.makePlot(sensSoln,'marginals',[],false,...
                    [fig1,fig2,fig3],{'b','linewidth',2})
% STL1 Model:
fig4 = figure(4);clf; set(fig4,'Name','Marginal Sensitivity, offGene');
fig5 = figure(5);clf; set(fig5,'Name','Marginal Sensitivity, onGene');
fig6 = figure(6);clf; set(fig6,'Name','Marginal Sensitivity, mRNA');
STL1_sens.makePlot(STL1_sensSoln,'marginals',[],false,...
                   [fig4,fig5,fig6],{'b','linewidth',2})

%% Compute FIMs using FSP sensitivity results
%% Model:
% Compute the FIM
Model_FIM = Model_sens;
fimResults = Model_FIM.computeFIM(sensSoln.sens); 

% Generate a count of measured cells (in place of real data)
cellCounts = 10*ones(size(Model_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[fimTotal,mleCovEstimate,fimMetrics] = ...
    Model_FIM.evaluateExperiment(fimResults,cellCounts)

% Plot the FIMs
fig7 = figure(7);clf; set(fig7,'Name',...
    'Fim-Predicted Uncertainty Ellipses');
Model_FIM.plotMHResults([],fimTotal,'log',[],fig7)
legend('FIM')

%% STL1 Model:
% Compute the FIM
STL1_FIM = STL1_sens;
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); 

% Generate a count of measured cells (i.e., in the case of real data,  
% the number of cells measured in the experiment)
STL1_cellCounts = 10*ones(size(STL1_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)

% Notice the difference in fimMetrics between the base model and STL1:

% fimMetrics = 
% 
%   struct with fields:
% 
%           det: 1.050283981892288e+12
%         trace: 1.875101485843730e+04
%     minEigVal: 1.010862869918453e+02

% STL1_fimMetrics = 
% 
%   struct with fields:
% 
%           det: 0
%         trace: 5.569716495312674e+04
%     minEigVal: 0

% This likely indicates a lack of identifiability between parameters in the
% STL1 Model, which unlike the base model has a time-varying input signal. 

% Let's try to tweak the STL1 time-varying input model a bit:

%% 
STL1 = SSIT;    

% Set species names for bursting gene expression model:
STL1.species = {'offGene';'onGene';'mRNA'}; 

% Set stoichiometry of reactions:
STL1.stoichiometry = [-1,0,0;...
                       1,0,0;...
                       0,1,-1]; 

% Define the time-varying TF/MAPK input signal:
STL1.inputExpressions = {'IHog',...
                              '(a0+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0))'};

% Set propensity functions:
STL1.propensityFunctions = {'IHog * offGene';'kr * onGene';'gr * mRNA'}; 

% Set initial guesses for parameters:
STL1.parameters = ({'kr',1; 'gr',0.5; ...
                    'a0',0.01; 'a1',1; 'r1',0.4; 'r2',.1});

% Set initial condition (one 'offGene'):
STL1.initialCondition = [1;0;0]; 


% Print a summary of STL1 Model:
STL1.summarizeModel


STL1.tSpan = linspace(0,20,200);

    % Create a copy of the STL1 Model for FSP:
    STL1_FSP = STL1;
    
    % Ensure the solution scheme is set to FSP (default):
    STL1_FSP.solutionScheme = 'FSP';  

    % This function compiles and stores the given reaction propensities  
    % into symbolic expression functions that use sparse matrices to  
    % operate on the system based on the current state. The functions are 
    % stored with the given prefix, in this case, 'STL1Model_FSP'
    STL1_FSP = STL1_FSP.formPropensitiesGeneral('STL1_FSP');
    
    % Set FSP 1-norm error tolerance:
    STL1_FSP.fspOptions.fspTol = 1e-4; 
    
    % Guess initial bounds on FSP StateSpace:
    STL1_FSP.fspOptions.bounds = [2,2,400];
    
    % Have FSP approximate the steady state for the initial distribution 
    % by finding the eigenvector corresponding to the smallest magnitude 
    % eigenvalue (i.e., zero, for generator matrix A, d/dtP(t)=AP(t)):
    STL1_FSP.fspOptions.initApproxSS = false; 
    
    % Solve Model:
    [STL1_FSPsoln,STL1_FSP.fspOptions.bounds] = STL1_FSP.solve; 
    
    % Plot marginal distributions:
    STL1_FSP.makePlot(STL1_FSPsoln,'marginals',[1:100:100],...
                       false,[1,2,3],{'linewidth',2})  
    STL1_FSP.makePlot(STL1_FSPsoln,'margmovie',[],false,[101],...
                       {'linewidth',2},'movie.mp4',[1,1,0.6],[2,3]) 

    STL1_FSP.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem
[STL1_sensSoln,STL1_bounds] = STL1_FSP.solve(STL1_FSPsoln.stateSpace); 

% STL1 Model:
fig4 = figure(4);clf; set(fig4,'Name','Marginal Sensitivity, offGene');
fig5 = figure(5);clf; set(fig5,'Name','Marginal Sensitivity, onGene');
fig6 = figure(6);clf; set(fig6,'Name','Marginal Sensitivity, mRNA');
STL1_FSP.makePlot(STL1_sensSoln,'marginals',[],false,...
                   [fig4,fig5,fig6],{'b','linewidth',2})

STL1_FIM = STL1_FSP;
STL1_fimResults = STL1_FIM.computeFIM(STL1_sensSoln.sens); 

% Generate a count of measured cells (i.e., in the case of real data,  
% the number of cells measured in the experiment)
STL1_cellCounts = 10*ones(size(STL1_FIM.tSpan));

% Evaluate the provided experiment design (in "cellCounts") 
% and produce an array of FIMs (one for each parameter set)
[STL1_fimTotal,STL1_mleCovEstimate,STL1_fimMetrics] = ...
    STL1_FIM.evaluateExperiment(STL1_fimResults,STL1_cellCounts)


% Plot the FIMs
fig8 = figure(5);clf; set(fig8,'Name',...
     'Fim-Predicted Uncertainty Ellipses');
STL1_FIM.plotMHResults([],STL1_fimTotal,'log',[],fig8)
legend('FIM')