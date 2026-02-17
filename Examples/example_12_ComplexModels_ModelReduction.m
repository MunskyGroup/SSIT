%% SSIT/Examples/example_12_ModelReduction
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
% addpath(genpath('../src'));

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
%   'Linear State Lumping' (LNSL) --
%       State space is divided into linearly distributed bins.
%       Then within each bin, the proability distribution is assumed to
%       be constant. The number of bins for each species is defined by
%       reductionOrder.
%   'Logarithmic State Lumping' (LGSL) --
%       State space is divided into logarithmically distributed bins.
%       Then within each bin, the proability distribution is assumed to
%       be constant. The number of bins for each species is defined by
%       reductionOrder.
%   'Principle Orthogonal Decomposition' (POD) - 
%       The infintesimal generator is projected onto the orthonorml
%       basis spanned by a previous solution of the full master
%       equation at discrete time points.  This is done by perfoming
%       SVD on the solutions and choosing the output space
%       corresponding to the 'redOrder' largest singular values.
%       For best use, this reduction should be found using a
%       fine time resolution in the calculation of the FSP. The size of the
%       reduced model must be specified as 'reductionOrder'. Because the
%       POD requires a solution of the FSP, this reduction is usually only
%       helpful for situations where many solutions are needed (e.g.,
%       during model fitting).
%   'No Transform' - 
%       Test case where no reduction is applied.
%   Other Options that are only aplicable to time invariant systems:
%   'Log Lump QSSA' (LGQSSA) -- 
%       Same as LGSL, but where the distribution within each lump is
%       assumed to be distributed according to the quasi-steady state
%       assumption. 
%   'Dynamic Mode Decomposition' (DM) -- 
%   (POD2) --
%       Extension of PDO that also includes the time derivatives of the CME
%       (i.e., A*P(t_i)) solution in the set of vectors onto which the CME
%       is projected. Only the vectors corresponding top the top
%       'reductionOrder' singular values are kept for the projection.
%   'Eigen Decomposition'(ED) --
%       The infintesimal generator is projected onto the eigenvectors
%       corresponding to its least negative eigenvalues. The number of
%       eigenvectors for this projection is defined by 'reductionOrder'.
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reductionType = 'POD'; 
reductionOrder = 50;
% podTimeSetSize = 30;    % Only needed for the POD reduction scheme.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Proper Orthogonal Decomposition (POD) to create a reduced 
%% model for computing FSP solutions for the 4-state time-varying 
%% STL1 yeast model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4-state STL1:
% Set up to solve FSP solution again to time following expansion:
STL1_MR_setup.fspOptions.initApproxSS = true;
STL1_MR_setup.tSpan = linspace(0,50,30);

% Print the computation time to solve the FSP using "tic" and "toc":
tic
% [STL1_FSPsoln_expand,STL1_MR_setup.fspOptions.bounds] = STL1_MR_setup.solve;
[~,~,STL1_MR_setup] = STL1_MR_setup.solve;
STL1_SolveTime = toc

% Turn off further FSP expansion:
STL1_MR_setup.fspOptions.fspTol = inf;

%% Solve the POD reduced model for STL1:
% Make a copy of the full model:
STL1_MR = STL1_MR_setup;

% Set the solver to use ModelReduction:
STL1_MR.modelReductionOptions.useModReduction = true;
% FSP expansion should be suppressed when using Model Reductions

% Set type and order of Model Reduction:
STL1_MR.modelReductionOptions.reductionType = reductionType;
STL1_MR.modelReductionOptions.reductionOrder = reductionOrder;

% Call the SSIT to compute the model reduction transformation matrices:
STL1_MR = STL1_MR.computeModelReductionTransformMatrices;

% Solve the reduced model:
tic
[~,~,STL1_MR] = STL1_MR.solve;
STL1_SolveTimeReduced = toc

%% Plot the full and reduced FSP solutions:
STL1_MR_setup.plotFSP([],STL1_MR_setup.species(5), 'means', [], [], {'linewidth',4},...
    Title='4-state STL1 (FSP full)', TitleFontSize=24, AxisLabelSize=18,...
    TickLabelSize=18, XLabel='Time', YLabel='Molecule Count',...
    LegendFontSize=15, LegendLocation='northeast', Colors=[0.23,0.67,0.2],...
    YLim = [0,32]);

STL1_MR.plotFSP([],STL1_MR.species(5), 'means', [], [], {'linewidth',4}, XLabel='Time',...
    YLabel='Molecule Count', Title='4-state STL1 (FSP reduced)',...
    TitleFontSize=24, AxisLabelSize=18, TickLabelSize=18,...
    Colors=[0.23,0.67,0.2], LegendFontSize=15, LegendLocation='northeast',...
    YLim = [0,32]);

%% Save reduced model & solution
% saveNames = unique({
%     'STL1_MR_setup'
%     'STL1_MR'
%     'STL1_FSPsoln_expand'
%     'STL1_FSPsoln_Red'
%     });
% 
% save('example_12_ComplexModels_ModelReduction',saveNames{:})