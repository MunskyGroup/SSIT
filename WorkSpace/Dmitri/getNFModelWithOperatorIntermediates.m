function model = getNFModelWithOperatorIntermediates(...
    M, initialCondition, parameters)
    arguments
        M (1,1) double {mustBeInteger, mustBePositive}
        initialCondition (3,1) double {mustBeInteger, mustBeNonnegative}
        parameters (5,2) cell
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
    
    % Add reactions:

    for dEmptyCntr = 1:M
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

    stoichDegradation = zeros(size(model.species, 1), 1);
    stoichDegradation(M + 1) = -1; % One X is consumed.
    model = model.addReaction('Kminus*X', stoichDegradation);

    stoichActivation = zeros(size(model.species, 1), 1);
    stoichActivation(1) = 1; % One Dempty1 is produced.
    stoichActivation(M + 2) = -1; % One Dfull is consumed.
    model = model.addReaction('Kon*Dfull', stoichActivation);

    model.initialCondition = initialCondition;
    model.parameters = parameters;
    
    model = model.formPropensitiesGeneral(['NFOperator', num2str(M)]);
end