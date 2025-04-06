%% example_8_ModelReduction
% Example script to demonstrate how to create reduced FSP models using 
% different types of projection-based transformations.
addpath('../../');

%% Preliminaries
% Load our models from example_1_CreateSSITModels; 
% Compute FSP solutions using example_2_SolveSSITModels_FSP
%% Comment out the following 3 lines if example_2_SolveSSITModels_FSP 
%% has already been run:
close all 
clear
example_2_SolveSSITModels_FSP

% View model summaries
Model_FSP.summarizeModel
STL1_FSP.summarizeModel

% Make new copies of our models
Model_MR_setup = Model_FSP;
STL1_MR_setup = STL1_FSP;

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

reductionType = 'Proper Orthogonal Decomposition'; 
reductionOrder = 100;
%qssaSpecies = 2;       % Only needed for the QSSA reduction scheme.
podTimeSetSize = 100;   % Only needed for the POD reduction scheme.

% Time varying model (DUSP1)
        % Model1 = SSIT;
        % Model1.species = {'ActiveGene';'mRNA'};
        % Model1.initialCondition = [0;0];
        % Model1.propensityFunctions = {'kon*(1+IGR)*(2-ActiveGene)';'koff*ActiveGene';'kr*ActiveGene';'gr*mRNA'};
        % Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
        % Model1.inputExpressions = {'IGR','a1*exp(-r1*t)*(1-exp(-r2*t))'};
        % Model1.parameters = ({'koff',0.14;'kon',0.14;'kr',10;'gr',0.01;...
        %     'a1',0.4;'r1',0.04;'r2',0.1});
        % Model1.fspOptions.initApproxSS = true;
        % Model1.tSpan = linspace(0,180,12);

%% Model:
% Solve FSP solution again to time following expansion 
Model_MR_setup.fspOptions.initApproxSS = true;
Model_MR_setup.tSpan = linspace(0,180,12);

tic
[Model_FSPsoln_expand,Model_MR_setup.fspOptions.bounds] = ...
    Model_MR_setup.solve(Model_FSPsoln.stateSpace);
Model_SolveTime = toc

% Turn off further FSP expansion
Model_MR_setup.fspOptions.fspTol = inf;

% If using the POD, we will also need to generate basis set using solution
% at finer resolution. Note -- this means that the POD will be inefficient
% for the initial set up of the reduction.  The benefits typically come
% from solving the model multiple times with different parameters sets.
if strcmp(reductionType,'Proper Orthogonal Decomposition')
    tSpan = Model_MR_setup.tSpan;
    Model_MR_setup.tSpan = linspace(min(Model_MR_setup.tSpan),...
        max(Model_MR_setup.tSpan),podTimeSetSize);
    [Model_FSPsoln_fine,Model_MR_setup.fspOptions.bounds] = ...
        Model_MR_setup.solve(Model_FSPsoln_expand.stateSpace);
    Model_MR_setup.tSpan = tSpan;
end

%% Solving the reduced models
% Make a copy of the full model
Model_MR = Model_MR_setup;

% Set the solver to use ModelReduction
Model_MR.modelReductionOptions.useModReduction = true;
% FSP expansion should be supressed when using Model Reductions

% Set type and order of Model Reduction
Model_MR.modelReductionOptions.reductionType = reductionType;
Model_MR.modelReductionOptions.reductionOrder = reductionOrder;
%Model_MR.modelReductionOptions.qssaSpecies = qssaSpecies;

% Call SSIT to compute the model reduction transformation matrices
Model_MR = ...
    Model_MR.computeModelReductionTransformMatrices(Model_FSPsoln_fine);

% Solve the reduce model
tic
Model_FSPsoln_Red = Model_MR.solve(Model_FSPsoln_fine.stateSpace);
Model_SolveTimeReduced = toc

% Make Figures to compare the results. Here, we will plot the original
% model in blue and the reduced model in red lines.
Model_MR_setup.makePlot(Model_FSPsoln_expand,'meansAndDevs',...
                        [],[],1,{'Color',[0,0,1]})
Model_MR.makePlot(Model_FSPsoln_Red,'meansAndDevs',...
                        [],[],1,{'Color',[1,0,0]})
figure(1);legend('Full','Reduced','Location','southeast')

Model_MR_setup.makePlot(Model_FSPsoln_expand,'marginals',...
                        [],[],[],{'Color',[0,0,1]})
Model_MR.makePlot(Model_FSPsoln_Red,'marginals',[],[],[],{'Color',[1,0,0]})
figure(2);legend('Full','Reduced','Location','eastoutside')
figure(3);legend('Full','Reduced','Location','eastoutside')

%% STL1 Model:

% Solve FSP solution again to time following expansion 
tic
[STL1_FSPsoln_expand,STL1_MR_setup.fspOptions.bounds] = ...
    STL1_MR_setup.solve(STL1_FSPsoln.stateSpace);
STL1_SolveTime = toc

% Turn off further FSP expansion
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

%% Solving the reduced models
% Make a copy of the full model.
STL1_MR = STL1_MR_setup;

% Set the solver to use ModelReduction
STL1_MR.modelReductionOptions.useModReduction = true;
% FSP expansion should be supressed when using Model Reductions

% Set type and order of Model Reduction
STL1_MR.modelReductionOptions.reductionType = reductionType;
STL1_MR.modelReductionOptions.reductionOrder = reductionOrder;
%Model_MR.modelReductionOptions.qssaSpecies = qssaSpecies;

% Call SSIT to compute the model reduction transformation matrices
STL1_MR = ...
    STL1_MR.computeModelReductionTransformMatrices(STL1_FSPsoln_fine);

% Solve the reduce model.
tic
STL1_FSPsoln_Red = STL1_MR.solve(STL1_FSPsoln_fine.stateSpace);
STL1_SolveTimeReduced = toc

% Make Figures to compare the results. Here, we will plot the original
% model in blue and the reduced model in red lines.
STL1_MR_setup.makePlot(STL1_FSPsoln_expand,'meansAndDevs',...
                        [],[],1,{'Color',[0,0,1]})
STL1_MR.makePlot(STL1_FSPsoln_Red,'meansAndDevs',...
                        [],[],1,{'Color',[1,0,0]})
figure(4);legend('Full','Reduced','Location','southeast')

STL1_MR_setup.makePlot(STL1_FSPsoln_expand,'marginals',...
                        [],[],[],{'Color',[0,0,1]})
STL1_MR.makePlot(STL1_FSPsoln_Red,'marginals',[],[],[],{'Color',[1,0,0]})
figure(5);legend('Full','Reduced','Location','eastoutside')
figure(6);legend('Full','Reduced','Location','eastoutside')