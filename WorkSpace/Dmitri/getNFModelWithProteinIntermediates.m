function model = getNFModelWithProteinIntermediates(M)
    arguments
        M (1,1) integer {mustBePositive}
    end
    
    model = SSIT;
    
    % Add all species, including the protein intermediates.
    % Each Xi will have species number i; X will therefore have number
    % M+1; Dempty, M+2; Dfull, M+3.

    for xCntr = 1:M
        model = model.addSpecies(['X', num2str(xCntr)], 0);
    end
    model = model.addSpecies('X', 0);
    model = model.addSpecies('Dempty', 1);
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

    for xCntr = 1:M-1
        stoichForwardProgress = zeros(size(model.species, 1), 1);
        % Each forward-progress reaction produces one Xi and consumes one
        % X(i+1).
        stoichForwardProgress([xCntr, xCntr + 1]) = [1, -1];
        model = model.addReaction(...
            ['M*X', num2str(xCntr + 1), '*', num2str(M), '/tau'], ...
            stoichForwardProgress);
    end

    stoichProduction = zeros(size(model.species, 1), 1);
    stoichProduction(M) = 1; % One XM is produced.
    model = model.addReaction('Kplus*Dempty', stoichProduction);

    stoichDegradation = zeros(size(model.species, 1), 1);
    stoichDegradation(M + 1) = -1; % One X is consumed.
    model = model.addReaction('Kminus*X', stoichDegradation);

    stoichActivation = zeros(size(model.species, 1), 1);
    stoichActivation([M + 2, M + 3]) = [1, -1]; % Dfull -> Dempty.
    model = model.addReaction('Kon*Dfull', stoichActivation);

    stoichDeactivation = zeros(size(model.species, 1), 1);
    stoichDeactivation([M + 2, M + 3]) = [-1, 1]; % Dempty -> Dfull.
    model = model.addReaction('Koff*Dempty*X', stoichDeactivation);

    model.summarizeModel
end