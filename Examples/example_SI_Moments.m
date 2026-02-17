%% SSIT/Examples/example_SI_Moments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2.2: Finding and visualizing master equation solutions
%   * Compute moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% clear
% close all
% addpath(genpath('../src'));

% example_1_CreateSSITModels
% example_2_SolveSSITModels_ODE

% Load the models created in example_1_CreateSSITModels:
% load('example_2_SolveSSITModels_ODE.mat')

% View model summaries:
Model_ODE.summarizeModel
STL1_ODE.summarizeModel
STL1_4state_ODE.summarizeModel

% Set the times at which distributions will be computed:
Model_ODE.tSpan = linspace(0,50,101);
STL1_ODE.tSpan = linspace(0,50,101);
STL1_4state_ODE.tSpan = linspace(0,50,101);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Bursting Gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the Bursting Gene model for moment equations
Model_mom = Model_ODE;
Model_mom.solutionScheme = 'moments';
[~,~,Model_mom] = Model_mom.solve;

% Number of species
nSp = numel(Model_mom.species);

%% First moments from the "moments" solver
% Moments: [ (#means + #secondMoments) x nTimes ]
mom = Model_mom.Solutions.moments;

% First nSp rows are the means ⟨x_i⟩
means_mom = mom(1:nSp, :);   % size: nSp x nTimes

%% Species trajectories from the ODE solver
% Get ODE solution:
Y_ode = Model_mom.Solutions.ode;   % nTimes x nSpecies
means_ode = Y_ode.';                     % transpose: nSp x nTimes

%% Compare ODEs vs moments 
% Compare mRNA:
i_mRNA = find(strcmp(Model_mom.species,'mRNA'));

den_mRNA = max(abs(means_ode(i_mRNA,:)), 1e-12);
errMean_mRNA = ...
    max(abs(means_mom(i_mRNA,:) - means_ode(i_mRNA,:)) ./ den_mRNA);

% Compare all species:
den_all = max(abs(means_ode), 1e-12);
errMean_all = max(max(abs(means_mom - means_ode) ./ den_all));

tol = 1e-2;   % 1% tolerance

if isequal(errMean_mRNA<tol,true)
    disp('Bursting Gene: mRNA mean from moments matches ODE within 1%');
else 
    disp('Bursting Gene: mRNA mean from moments does not match ODE within 1%');
end

if isequal(errMean_all<tol,true)
    disp('Bursting Gene: species means from moments match ODE within 1%');
else 
    disp('Bursting Gene: species means from moments do not match ODE within 1%');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): STL1 (simple)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the (simple) STL1 model for moment equations
STL1_mom = STL1_ODE;
STL1_mom.solutionScheme = 'moments';
[~,~,STL1_mom] = STL1_mom.solve;

% Number of species
nSp = numel(STL1_mom.species);

%% First moments from the "moments" solver
% Moments: [ (#means + #secondMoments) x nTimes ]
mom = STL1_mom.Solutions.moments;

% First nSp rows are the means ⟨x_i⟩
means_mom = mom(1:nSp, :);   % size: nSp x nTimes

%% Species trajectories from the ODE solver
% Get ODE solution:
Y_ode = STL1_mom.Solutions.ode;   % nTimes x nSpecies
means_ode = Y_ode.';                     % transpose: nSp x nTimes

%% Compare ODEs vs moments 
% Compare mRNA:
i_mRNA = find(strcmp(STL1_mom.species,'mRNA'));

den_mRNA = max(abs(means_ode(i_mRNA,:)), 1e-12);
errMean_mRNA = ...
    max(abs(means_mom(i_mRNA,:) - means_ode(i_mRNA,:)) ./ den_mRNA);

% Compare all species:
den_all = max(abs(means_ode), 1e-12);
errMean_all = max(max(abs(means_mom - means_ode) ./ den_all));

tol = 1e-2;   % 1% tolerance

if isequal(errMean_mRNA<tol,true)
    disp('STL1: mRNA mean from moments matches ODE within 1%');
else 
    disp('STL1: mRNA mean from moments does not match ODE within 1%');
end

if isequal(errMean_all<tol,true)
    disp('STL1: species means from moments match ODE within 1%');
else 
    disp('STL1: species means from moments do not match ODE within 1%');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): 4-state STL1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a copy of the 4-state STL1 model for moment equations
STL1_4state_mom = STL1_4state_ODE;
STL1_4state_mom.solutionScheme = 'moments';
[~,~,STL1_4state_mom] = STL1_4state_mom.solve;

% Number of species
nSp = numel(STL1_4state_mom.species);

%% First moments from the "moments" solver
% Moments: [ (#means + #secondMoments) x nTimes ]
mom = STL1_4state_mom.Solutions.moments;

% First nSp rows are the means ⟨x_i⟩
means_mom = mom(1:nSp, :);   % size: nSp x nTimes

%% Species trajectories from the ODE solver
% Get ODE solution:
Y_ode = STL1_4state_ODE.Solutions.ode;   % nTimes x nSpecies
means_ode = Y_ode.';                     % transpose: nSp x nTimes

%% Compare ODEs vs moments 
% Compare mRNA:
i_mRNA = find(strcmp(STL1_4state_mom.species,'mRNA'));

den_mRNA = max(abs(means_ode(i_mRNA,:)), 1e-12);
errMean_mRNA = ...
    max(abs(means_mom(i_mRNA,:) - means_ode(i_mRNA,:)) ./ den_mRNA);

% Compare all species:
den_all = max(abs(means_ode), 1e-12);
errMean_all = max(max(abs(means_mom - means_ode) ./ den_all));

tol = 1e-2;   % 1% tolerance

if isequal(errMean_mRNA<tol,true)
    disp('4-state STL1: mRNA mean from moments matches ODE within 1%');
else 
    disp('4-state STL1: mRNA mean from moments does not match ODE within 1%');
end

if isequal(errMean_all<tol,true)
    disp('4-state STL1: species means from moments match ODE within 1%');
else 
    disp('4-state STL1: species means from moments do not match ODE within 1%');
end


