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

initialConditionX = 44; % Result from actual SSA runs
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

datasetSizes = [1, 2, 5, 10, 15, 20, 30] * 1e3;
datasetFractions = datasetSizes / 3e4;
for dsCntr = 1:length(datasetSizes)
    datasetNames{dsCntr} = ['ssa_', num2str(datasetSizes(dsCntr)), '.csv'];
end

fspTimes = zeros(length(M), 1);
fitTimes = zeros(length(M), length(datasetSizes));

%% %% Protein intermediates

paramsToFit = [1, 4, 5]; % kPlus, kOff, tau

fitParamValues = zeros(length(M), length(datasetSizes), length(paramsToFit));

for Mcntr = 1:length(M)
    curM = M(Mcntr);
    piModel = getNFModelWithProteinIntermediates(...
        curM, initialConditionX, parameters);
    piModel.tSpan = tSpan;
    piModel.solutionScheme = 'FSP';
    tic;
    [piFSPsoln, piModel.fspOptions.bounds] = piModel.solve;
    fspTimes(Mcntr) = toc;

    for dsCntr = 1:length(datasetSizes)
        curModel = piModel; % Clone the model solved via FSP
        
        % Load the subset of the SSA data
        
        curModel = curModel.loadData(datasetNames{dsCntr}, ...
            {'Dempty','exp1_s1','Dfull','exp1_s2','X','exp1_s3'});

        % Only a subset of the parameters will be made "free"

        curModel.fittingOptions.modelVarsToFit = paramsToFit;

        fitOptions = optimset('Display', 'none', 'MaxIter', 1000);
        tic;
        fitParameters = curModel.maximizeLikelihood([], fitOptions);        
        fitTimes(Mcntr,dsCntr) = toc; 
        curModel.parameters(curModel.fittingOptions.modelVarsToFit, 2) = ...
            num2cell(fitParameters);
        for paramCntr = 1:length(paramsToFit)
            fitParamValues(Mcntr, dsCntr, paramCntr) = ...
                curModel.parameters{paramsToFit(paramCntr), 2};
        end
    end % for dsCntr = 1:length(datasetSizes) ...

    %piModel.makePlot(piFSPsoln,'meansAndDevs',[],[]) % Make plot of mean vs. time.
    %piModel.makePlot(piFSPsoln,'marginals',[],[]) % Make plot of marginals
end % for Mcntr = 1:length(M)

disp(fspTimes)
save piFSPtimes.mat fspTimes
disp(fitTimes)
save piFitTimes.mat fitTimes
disp(fitParamValues)
save piFitParamValues.mat fitParamValues