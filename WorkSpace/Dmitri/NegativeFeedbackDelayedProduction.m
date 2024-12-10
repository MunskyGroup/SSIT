%% Initial Setup

close all
clear all
addpath(genpath('../../src'));

kOn = 0.1;
kOff = kOn * 1000 / 2;

tau = 50;
delayPeriods = 30;
tSpan = linspace(0, tau * delayPeriods, 1 + delayPeriods);

Model = SSIT;
Model.species = {'Dempty'; 'Dfull'; 'X'};

initialCondition = [1; 0; randi(50)];
Model.initialCondition = initialCondition;

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
tic;
SSAsoln = Model.solve;
SSAtime = toc;

%% Plot trajectories
Model.makePlot(SSAsoln,'trajectories',[],[],5) % Make some plots.

%% Protein intermediates - 3

pi3 = getNFModelWithProteinIntermediates(3, initialCondition, parameters);

%% FSP Solution
pi3.tSpan = tSpan;
pi3.solutionScheme = 'FSP';
tic;
[pi3FSPsoln, pi3.fspOptions.bounds] = pi3.solve;
pi3FSPtime = toc;

%% Plot FSP solution
pi3.makePlot(pi3FSPsoln,'meansAndDevs',[],[],1,{'linewidth',3,'color',[0,1,1]}) % Make plot of mean vs. time.
pi3.makePlot(pi3FSPsoln,'marginals',[],[],2,{'linewidth',3,'color',[0,0,1]}) % Make plot of mean vs. time.

%% Operator intermediates - 3

oi3 = getNFModelWithOperatorIntermediates(3, initialCondition, parameters);

%% FSP solution
oi3.tSpan = tSpan;
oi3.solutionScheme = 'FSP';
tic;
[oi3FSPsoln, oi3.fspOptions.bounds] = oi3.solve;
oi3FSPtime = toc;

%% Plot FSP solution
oi3.makePlot(oi3FSPsoln,'meansAndDevs',[],[],3,{'linewidth',3,'color',[0,1,1]}) % Make plot of mean vs. time.
oi3.makePlot(oi3FSPsoln,'marginals',[],[],4,{'linewidth',3,'color',[0,0,1]}) % Make plot of mean vs. time.
