function model = getNFModelWithProteinIntermediates(...
    M, initialConditionX, parameters)
    arguments
        M (1,1) double {mustBeInteger, mustBePositive}
        initialConditionX (1,1) double {mustBeInteger, mustBeNonnegative}
        parameters (5,2) cell
    end
    
    model = SSIT;
    model.species = {};
    model.initialCondition = [];
    model.propensityFunctions = {};
    model.stoichiometry = [];
    model.parameters = {};
    
    % Add all species, including the protein intermediates.
    % Each Xi will have species number i; X will therefore have number
    % M+1; Dempty, M+2; Dfull, M+3.

    for xCntr = 1:M
        model = model.addSpecies(['X', num2str(xCntr)], 0);
    end
    model = model.addSpecies('X', initialConditionX);
    model = model.addSpecies('Dempty', 1);
    model = model.addSpecies('Dfull', 0);

    model.stoichiometry = [];
    
    % Add reactions:

    for xCntr = 1:M-1
        stoichForwardProgress = zeros(size(model.species, 1), 1);
        % Each forward-progress reaction produces one Xi and consumes one
        % X(i+1).
        stoichForwardProgress([xCntr, xCntr + 1]) = [1, -1];
        model = model.addReaction(...
            ['X', num2str(xCntr + 1), '*', num2str(M), '/tau'], ...
            stoichForwardProgress);
    end

    stoichForwardProgress = zeros(size(model.species, 1), 1);
    % The final forward-progress reaction produces one X and consumes one
    % X1.
    stoichForwardProgress(1) = -1;
    stoichForwardProgress(M + 1) = 1;
    model = model.addReaction(['X1*', num2str(M), '/tau'], ...
        stoichForwardProgress);

    stoichProduction = zeros(size(model.species, 1), 1);
    stoichProduction(M) = 1; % One XM is produced.
    model = model.addReaction('kPlus*Dempty', stoichProduction);

    stoichDegradation = zeros(size(model.species, 1), 1);
    stoichDegradation(M + 1) = -1; % One X is consumed.
    model = model.addReaction('kMinus*X', stoichDegradation);

    stoichActivation = zeros(size(model.species, 1), 1);
    stoichActivation([M + 2, M + 3]) = [1, -1]; % Dfull -> Dempty.
    model = model.addReaction('kOn*Dfull', stoichActivation);

    stoichDeactivation = zeros(size(model.species, 1), 1);
    stoichDeactivation([M + 2, M + 3]) = [-1, 1]; % Dempty -> Dfull.
    model = model.addReaction('kOff*Dempty*X', stoichDeactivation);

    model.parameters = parameters;
    
    model = model.formPropensitiesGeneral(['NFProtein', num2str(M)]);
end