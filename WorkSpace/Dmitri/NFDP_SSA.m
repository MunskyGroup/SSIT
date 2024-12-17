%% Initial Setup

close all
clear all
addpath(genpath('../../src'));

M = [1 2 3 6]; % Number of intermediate states

kOn = 0.1;
kOff = kOn * 1000 / 2;

tau = 50;
delayPeriods = 10;
tSpan = linspace(0, tau * delayPeriods, 1 + delayPeriods);

Model = SSIT;
Model.species = {'Dempty'; 'Dfull'; 'X'};

initialConditionX = randi([15, 50]);
Model.initialCondition = [1; 0; initialConditionX];

parameters = ({...
    'kPlus', 100; ... % (Delayed) rate of protein production
    'kMinus', 1; ... % Rate of protein degradation
    'kOn', kOn; ... % Rate of operator site emptying / activation
    'kOff', kOff; ... % Rate of operator site filling / deactivation
    'tau', tau ... % Delay in protein production
    });
Model.parameters = parameters;

Model.propensityFunctions = {
    'kOn*Dfull', ... % Operator site emptying
    'kOff*Dempty*X', ... % Operator site filling
    'kPlus*Dempty', ... % (Delayed) protein production
    'kMinus*X' ... % Protein degradation
    };

Model.stoichiometry = [...
     1, -1,  0,  0; ... % Dempty
    -1,  1,  0,  0; ... % Dfull
     0,  0,  1, -1 ... % X
    ];

% The only delayed reaction in the system is the third, with a delay equal
% to tau. It has a single reactant, one Dempty, so we must include it in
% the delayed-reaction scheduling stochiometries:

Model.delayedReactions = [3, Model.parameters{5,2}];

Model.delayedReactionSchedulingS = [...
     0,  0, -1,  0; ... % Dempty
     0,  0,  0,  0; ... % Dfull
     0,  0,  0,  0 ... % X
    ];

%% Configure and run SSA

Model.tSpan = tSpan;
Model = Model.formPropensitiesGeneral('NFDP');

Model.solutionScheme = 'SSA';

% Number of independent data sets to generate.
Model.ssaOptions.Nexp = 2;   

% Number of cells to include at each time point for each data set.
Model.ssaOptions.nSimsPerExpt = 3e4; 

tic;
SSAsoln = Model.solve([], 'ssa_out.csv');
SSAtime = toc;

save SSAtime.mat SSAtime
save SSAsoln.mat SSAsoln