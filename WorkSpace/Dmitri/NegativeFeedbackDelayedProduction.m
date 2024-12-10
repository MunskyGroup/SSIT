%% Initial Setup

close all
clear all
addpath(genpath('../../src'));

M = [1 2 3 6 10]; % Number of intermediate states

kOn = 0.1;
kOff = kOn * 1000 / 2;

tau = 50;
delayPeriods = 10;
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
Model.ssaOptions.Nexp = 2;   % Number of independent data sets to generate.
Model.ssaOptions.nSimsPerExpt = 1000; % Number of cells to include at each time point for each data set.
tic;
SSAsoln = Model.solve;
SSAtime = toc;

save("SSAsoln", SSAsoln);

%% Plot trajectories
Model.makePlot(SSAsoln,'trajectories',[],[]) % Make some plots.

%% Protein intermediates

piFSPtimes = zerosLike(M);
for Mcntr = 1:length(M)
    curM = M(Mcntr);
    piModel = getNFModelWithProteinIntermediates(...
        curM, initialCondition, parameters);
    piModel.tSpan = tSpan;
    piModel.solutionScheme = 'FSP';
    tic;
    [piFSPsoln, piModel.fspOptions.bounds] = piModel.solve;
    piFSPtimes(Mcntr) = toc;

    piModel.makePlot(piFSPsoln,'meansAndDevs',[],[]) % Make plot of mean vs. time.
    piModel.makePlot(piFSPsoln,'marginals',[],[]) % Make plot of marginals
end

%% Operator intermediates

opFSPtimes = zerosLike(M);
for Mcntr = 1:length(M)
    curM = M(Mcntr);
    opModel = getNFModelWithOperatorIntermediates(...
        curM, initialCondition, parameters);
    opModel.tSpan = tSpan;
    opModel.solutionScheme = 'FSP';
    tic;
    [piFSPsoln, opModel.fspOptions.bounds] = opModel.solve;
    opFSPtimes(Mcntr) = toc;

    opModel.makePlot(piFSPsoln,'meansAndDevs',[],[]) % Make plot of mean vs. time.
    opModel.makePlot(piFSPsoln,'marginals',[],[]) % Make plot of marginals
end

disp(piFSPtimes)
save("piFSPtimes", piFSPtimes);
disp(opFSPtimes)
save("opFSPtimes", opFSPtimes);