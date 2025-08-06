%% example_6_SensitivityAnalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.3: Sensitivity analysis and Fisher Information Matrix (Part I)
%   * Compute model sensitivities to small changes in the model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_1_CreateSSITModels and  
% computed FSP solutions from example_4_SolveSSITModels_FSP

% clear
% close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels  
% example_4_SolveSSITModels_FSP

loadPrevious = true;
savedWorkspace = 'example_6_SensitivityAnalysis';

% View model summaries
Model_FSP.summarizeModel
STL1_FSP.summarizeModel
STL1_4state_FSP.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Solve sensitivities of the bursting gene model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the bursting gene model solved by FSP for sensitivity 
% analysis:
Model_sens = Model_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
Model_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem:
[Model_sensSoln,Model_bounds] = Model_sens.solve(Model_FSPsoln.stateSpace);

% Plot the results from the sensitivity analysis:
fig1 = figure(1);clf; set(fig1,'Name','Marginal Sensitivity, offGene');
fig2 = figure(2);clf; set(fig2,'Name','Marginal Sensitivity, onGene');
fig3 = figure(3);clf; set(fig3,'Name','Marginal Sensitivity, mRNA');
Model_sens.makePlot(Model_sensSoln,'marginals',[],false,...
                    [fig1,fig2,fig3],{'b','linewidth',2})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve sensitivities of the time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the time-varying STL1 yeast model solved by FSP for 
% sensitivity analysis:
STL1_sens = STL1_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
STL1_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem: 
[STL1_sensSoln,STL1_bounds] = STL1_sens.solve(STL1_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis:
fig4 = figure(4);clf; set(fig4,'Name','Marginal Sensitivity, offGene');
fig5 = figure(5);clf; set(fig5,'Name','Marginal Sensitivity, onGene');
fig6 = figure(6);clf; set(fig6,'Name','Marginal Sensitivity, mRNA');
STL1_sens.makePlot(STL1_sensSoln,'marginals',[],false,...
                   [fig4,fig5,fig6],{'b','linewidth',2})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve sensitivities of the 4-state time-varying STL1 yeast model
%  from example_1_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the time-varying STL1 yeast model solved by FSP for 
% sensitivity analysis:
STL1_4state_sens = STL1_4state_FSP;

%% Solve FSP sensitivities
% Set solution schemes to FSP sensitivity:
STL1_4state_sens.solutionScheme = 'fspSens'; 

% Solve the sensitivity problem: 
[STL1_4state_sensSoln,STL1_4state_bounds] = ...
    STL1_4state_sens.solve(STL1_4state_FSPsoln.stateSpace); 

% Plot the results from the sensitivity analysis:
fig7 = figure(7);clf; set(fig7,'Name','Marginal Sensitivity, s1');
fig8 = figure(8);clf; set(fig8,'Name','Marginal Sensitivity, s2');
fig9 = figure(9);clf; set(fig9,'Name','Marginal Sensitivity, s3');
fig10 = figure(10);clf; set(fig10,'Name','Marginal Sensitivity, s4');
fig11 = figure(11);clf; set(fig11,'Name','Marginal Sensitivity, mRNA');
STL1_4state_sens.makePlot(STL1_4state_sensSoln,'marginals',[],false,...
                   [fig7,fig8,fig9,fig10,fig11],{'b','linewidth',2})

saveNames = unique({'Model_sens'
    'Model_sensSoln'
    'Model_bounds'
    'STL1_sens'
    'STL1_sensSoln'
    'STL1_bounds'
    'STL1_4state_sens'
    'STL1_4state_sensSoln'
    'STL1_4state_bounds'
    });
    
save('example_6_SensitivityAnalysis',saveNames{:})