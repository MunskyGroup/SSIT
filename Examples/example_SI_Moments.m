%% SSIT/Examples/example_SI_Moments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moment solutions
%   * Compute moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% clear
% close all
% addpath(genpath('../src'));

% example_01_CreateSSITModels
% example_02_SolveSSITModels_ODE

% Load the models created in example_01_CreateSSITModels:
load('ExampleSaveFiles/example_2_SolveSSITModels_ODE.mat')

% View model summaries:
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

% Set the times at which distributions will be computed:
Model.tSpan = linspace(0,50,101);
STL1.tSpan = linspace(0,50,101);
STL1_4state.tSpan = linspace(0,50,101);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Bursting Gene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model = Model.solve(solver='moments');

% Number of species
nSp = numel(Model.species);

%% First moments from the "moments" solver
% Moments: [ (#means + #secondMoments) x nTimes ]
mom = Model.Solutions.moments;

% First nSp rows are the means ⟨x_i⟩
means_mom = mom(1:nSp, :);   % size: nSp x nTimes

%% Species trajectories from the ODE solver
% Get ODE solution:
Y_ode = Model.Solutions.ode;   % nTimes x nSpecies
means_ode = Y_ode.';                     % transpose: nSp x nTimes

%% Compare ODEs vs moments 
% Compare mRNA:
i_mRNA = find(strcmp(Model.species,'mRNA'));

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

STL1 = STL1.solve(solver='moments');

% Number of species
nSp = numel(STL1.species);

%% First moments from the "moments" solver
% Moments: [ (#means + #secondMoments) x nTimes ]
mom = STL1.Solutions.moments;

% First nSp rows are the means ⟨x_i⟩
means_mom = mom(1:nSp, :);   % size: nSp x nTimes

%% Species trajectories from the ODE solver
% Get ODE solution:
Y_ode = STL1.Solutions.ode;   % nTimes x nSpecies
means_ode = Y_ode.';                     % transpose: nSp x nTimes

%% Compare ODEs vs moments 
% Compare mRNA:
i_mRNA = find(strcmp(STL1.species,'mRNA'));

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
STL1_4state = STL1_4state.solve(solver='moments');

% Number of species
nSp = numel(STL1_4state.species);

%% First moments from the "moments" solver
% Moments: [ (#means + #secondMoments) x nTimes ]
mom = STL1_4state.Solutions.moments;

% First nSp rows are the means ⟨x_i⟩
means_mom = mom(1:nSp, :);   % size: nSp x nTimes

%% Species trajectories from the ODE solver
% Get ODE solution:
Y_ode = STL1_4state.Solutions.ode;   % nTimes x nSpecies
means_ode = Y_ode.';                     % transpose: nSp x nTimes

%% Compare ODEs vs moments 
% Compare mRNA:
i_mRNA = find(strcmp(STL1_4state.species,'mRNA'));

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

%% Plot results:
STL1_4state.plotMoments(solution=STL1_4state.Solutions.moments,...
    speciesNames=STL1_4state.species(5), plotType="meansanddevs",...
    indTimes=STL1_4state.tSpan, lineProps={'linewidth',4},...
    Title='4-state STL1 (mRNA)', TitleFontSize=26, LegendFontSize=20,...
    LegendLocation='east', Colors=[0.23,0.67,0.20],...
    YLabel='Molecule Count', AxisLabelSize=20, TickLabelSize=20);

STL1_4state.plotMoments(solution=STL1_4state.Solutions.moments,...
    speciesNames=STL1_4state.species(1:4), plotType="meansanddevs",...
    indTimes=STL1_4state.tSpan, lineProps={'linewidth',4},...
    Title='4-state STL1 (mRNA)', TitleFontSize=26,...
    LegendLocation='east', YLabel='Molecule Count',...
    AxisLabelSize=20, TickLabelSize=20, LegendFontSize=20);