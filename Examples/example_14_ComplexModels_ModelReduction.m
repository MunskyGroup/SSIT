%% example_14_ModelReduction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.5: Complex models
%   * Create reduced FSP models using different types of 
%     projection-based transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the STL1 model from example_1_CreateSSITModels and FSP solutions from 
% example_4_SolveSSITModels_FSP
% clear
% close all
addpath(genpath('../../src'));

% example_1_CreateSSITModels 
% example_4_SolveSSITModels_FSP
% example_8b_LoadingandFittingData_SimulatingData

%% Load pre-computed FSP solutions:
% load('example_4_SolveSSITModels_FSP.mat')

% View model summaries:
STL1_4state_FSP.summarizeModel

% Make a copy of the STL1 model to set up for model reduction:
STL1_MR_setup = STL1_4state_FSP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose which type of model reduction to apply. Options include:
%   'Proper Orthogonal Decomposition' - solve the FSP once and then use POD
%       to construct a reduced basis set that covers the current FSP
%       solution.  For best use, this reduction should be found using a
%       fine time resolution in the calculation of the FSP. The size of the
%       reduced model must be specified as 'reductionOrder'. Because the
%       POD requires a solution of the FSP, this reduction is usually only
%       helpful for situations where many solutions are needed (e.g.,
%       during model fitting).
%   'Log Lump QSSA' - forms a coarse rectangular mesh with grid points 
%       chosen logarithmically using the current FSP bounds. The number of 
%       grid lines must be specified in 'reductionOrder'.
%   'Eigen Decomposition Initial' - reduction to consider only the space
%       spanned by the initial condition plus the eigvenvectors 
%       corresponding to the eigenvalues with the largest real values.  
%       The number of modes to consider in the reduction is specified in 
%       'reductionOrder'. For time-varying systems, the basis vectors are 
%       found using the infinitesimal generator at t=0.  
%   'QSSA' - Reduction using QSSA applied to a specific species or set of
%       species. The list of species to be assumed at QSSA must be
%       specified in a vector 'reductionSpecies'. 
%   'No Transform' - test case where no reduction is applied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reductionType = 'Proper Orthogonal Decomposition'; 
reductionOrder = 100;
%qssaSpecies = 2;       % Only needed for the QSSA reduction scheme.
podTimeSetSize = 100;   % Only needed for the POD reduction scheme.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Proper Orthogonal Decomposition (POD) to create a reduced 
%% model for computing FSP solutions for the 4-state time-varying 
%% STL1 yeast model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STL1:
% Set up to solve FSP solution again to time following expansion:
STL1_MR_setup.fspOptions.initApproxSS = true;
STL1_MR_setup.tSpan = linspace(0,180,12);

% Print the computation time to solve the FSP using "tic" and "toc":
tic
[STL1_FSPsoln_expand,STL1_MR_setup.fspOptions.bounds] = ...
    STL1_MR_setup.solve(STL1_4state_FSPsoln.stateSpace);
STL1_SolveTime = toc

% Turn off further FSP expansion:
STL1_MR_setup.fspOptions.fspTol = inf;

% If using the POD, we will also need to generate basis set using solution
% at finer resolution. Note -- this means that the POD will be inefficient
% for the initial set up of the reduction.  The benefits typically come
% from solving the model multiple times with different parameters sets.
if strcmp(reductionType,'Proper Orthogonal Decomposition')
    tSpan = STL1_MR_setup.tSpan;
    STL1_MR_setup.tSpan = linspace(min(STL1_MR_setup.tSpan),...
        max(STL1_MR_setup.tSpan),podTimeSetSize);
    [STL1_FSPsoln_fine,STL1_MR_setup.fspOptions.bounds] = ...
        STL1_MR_setup.solve(STL1_FSPsoln_expand.stateSpace);
    STL1_MR_setup.tSpan = tSpan;
end

%% Solve the POD reduced model for STL1:
% Make a copy of the full model:
STL1_MR = STL1_MR_setup;

% Set the solver to use ModelReduction:
STL1_MR.modelReductionOptions.useModReduction = true;
% FSP expansion should be supressed when using Model Reductions

% Set type and order of Model Reduction:
STL1_MR.modelReductionOptions.reductionType = reductionType;
STL1_MR.modelReductionOptions.reductionOrder = reductionOrder;
%Model_MR.modelReductionOptions.qssaSpecies = qssaSpecies;

% Call SSIT to compute the model reduction transformation matrices:
STL1_MR = ...
    STL1_MR.computeModelReductionTransformMatrices(STL1_FSPsoln_fine);

% Solve the reduce model:
tic
STL1_FSPsoln_Red = STL1_MR.solve(STL1_FSPsoln_fine.stateSpace);
STL1_SolveTimeReduced = toc

% Make Figures to compare the results. Here, we will plot the original
% model in blue and the reduced model in red lines.
STL1_MR_setup.makePlot(STL1_FSPsoln_expand,'meansAndDevs',...
                        [],[],1,{'Color',[0,0,1]})
STL1_MR.makePlot(STL1_FSPsoln_Red,'meansAndDevs',...
                        [],[],1,{'Color',[1,0,0]})
figure(1);legend('Full','Reduced','Location','southeast')

STL1_MR_setup.makePlot(STL1_FSPsoln_expand,'marginals',...
                        [],[],[],{'Color',[0,0,1]})
STL1_MR.makePlot(STL1_FSPsoln_Red,'marginals',[],[],[],{'Color',[1,0,0]})
figure(2);legend('Full','Reduced','Location','eastoutside')
figure(3);legend('Full','Reduced','Location','eastoutside')