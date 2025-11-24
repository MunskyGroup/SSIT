%% SSIT/Examples/example_18_Moments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 
%   * 
%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% clear
% close all
% addpath(genpath('../src'));

% example_1_CreateSSITModels

% Load the models created in example_1_CreateSSITModels:
% load('example_1_CreateSSITModels.mat')

% View model summaries:
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

% Set the times at which distributions will be computed:
Model_FSP.tSpan = linspace(0,50,200);
STL1_FSP.tSpan = linspace(0,50,200);
STL1_4state_FSP.tSpan = linspace(0,50,200);


% Make a copy of the STL1 model for moment equations
STL1_4state_mom = STL1_4state;
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
    disp('mRNA mean from moments matches ODE within 1%');
else 
    disp('mRNA mean from moments does not match ODE within 1%');
end

if isequal(errMean_all<tol,true)
    disp('Means from moments matches ODE within 1%');
else 
    disp('Means from moments do not match ODE within 1%');
end


