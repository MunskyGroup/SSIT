%% Initial Setup

close all
clear all
addpath(genpath('../../src'));

Model = SSIT;
Model.species = {'Dempty'; 'Dfull'; 'X'};
Model.initialCondition = [1; 0; 0];
Model.parameters = ({...
    'kPlus', 100; ... % (Delayed) rate of protein production
    'kMinus', 1; ... % Rate of protein degradation
    'kMinus1', 2; ... % Rate of operator site emptying
    'k1', 1000; ... % Rate of operator site filling
    'tau', 2 ... % Delay in protein production
    });

Model.propensityFunctions = {
    'kMinus1*Dfull', ... % Operator site emptying
    'k1*Dempty*X', ... % Operator site filling
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

%Model.tSpan = linspace(0, 400, 401);
Model.tSpan = linspace(0, 9, 10);
Model = Model.formPropensitiesGeneral('NFDP');

Model.solutionScheme = 'SSA';
SSASoln = Model.solve;

%% Plot trajectories
Model.makePlot(SSASoln,'trajectories',[],[],4) % Make some plots.