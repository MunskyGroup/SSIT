function model = getNFModelWithOperatorIntermediates(M)
    arguments
        M (1,1) double {mustBeInteger, mustBePositive}
    end
    
    model = SSIT;
    model.species = {};
    model.initialCondition = [];
    model.propensityFunctions = {};
    model.stoichiometry = [];
    model.parameters = {};
    
    % Add all species, including the operator intermediates.
    % Each DemptyI will have species number I; X will therefore have number
    % M+1; Dfull, M+2.

    for dEmptyCntr = 1:M
        initialConditions = zeros(M, 1);
        initialConditions(1) = 1; % Only Dempty1 will be present initially.
        model = model.addSpecies(...
            ['Dempty', num2str(dEmptyCntr)], ...
            initialConditions(dEmptyCntr));
    end
    model = model.addSpecies('X', 0);    
    model = model.addSpecies('Dfull', 0);

    model.stoichiometry = [];

    % Add parameters:

    model.parameters = ({...
        'kPlus', 100; ... % (Delayed) rate of protein production
        'kMinus', 1; ... % Rate of protein degradation
        'kOn', 2; ... % Rate of operator site emptying
        'kOff', 1000; ... % Rate of operator site filling
        'tau', 50 ... % Delay in protein production
    });
    
    % Add reactions:

    for dEmptyCntr = 1:M-1
        % Each production reaction produces one X and shifts the operator
        % state, consuming one DemptyI and producing one Dempty(I+1).
        stoichProduction = zeros(size(model.species, 1), 1);
        stoichProduction([dEmptyCntr, dEmptyCntr + 1]) = [-1, 1];
        stoichProduction(M + 1) = 1;
        model = model.addReaction(...
            ['kPlus*Dempty', num2str(dEmptyCntr), ...
                '*', num2str(M), '/tau'], ...
            stoichProduction);
    end

    % The final production reaction produces one X and shifts the operator
    % state, consuming one DemptyM and producing one Dfull.

    stoichProduction = zeros(size(model.species, 1), 1);
    stoichProduction(M) = -1; % DemptyM
    stoichProduction(M + 1) = 1; % X
    stoichProduction(M + 2) = 1; % Dfull

    model = model.addReaction(...
        ['kPlus*Dempty', num2str(M), '*', num2str(M), '/tau'], ...
        stoichProduction);

    stoichDegradation = zeros(size(model.species, 1), 1);
    stoichDegradation(M + 1) = -1; % One X is consumed.
    model = model.addReaction('kMinus*X', stoichDegradation);

    stoichActivation = zeros(size(model.species, 1), 1);
    stoichActivation(1) = 1; % One Dempty1 is produced.
    stoichActivation(M + 2) = -1; % One Dfull is consumed.
    model = model.addReaction('kOn*Dfull', stoichActivation);

    model = model.formPropensitiesGeneral(['NFOperator', num2str(M)]);
end