function [X_array] = runSingleSsa(x0, S, W, T_array, isTimeVarying, ...
    signalUpdateRate, parameters, ...
    delayedReactions, delayedReactionSchedulingS)

% Start the simulation.

%% Step 1: Initialization
% Input values for the initial state from the provided vector, set time t
% to be the initial time of simulation, and set the reaction counter i = 1.

t = min(T_array);   % initial time of simulation.
x = x0;
iprint = 1;  % The next time at which to record results.
Nsp = size(S,1);  %Number of species.
Nt = length(T_array);
X_array = zeros(Nsp,Nt);
props = [];

hasDelayedReactions = ~isempty(delayedReactions);

% To avoid hints/warnings, we will initialize an array that will hold a
% record of delayed reactions: the first column will contain the scheduled
% time; the second, the reaction number.
scheduledDelayedReactions = [NaN, NaN];

if isTimeVarying
        S = [zeros(Nsp,1),S];
        props(1) = signalUpdateRate;
        jt = 1;
else
    jt=0;
end

isAFun = isa(W{1}, 'function_handle');
    
while t < max(T_array)
    %% Step 2: Compute propensities for all reactions
    if ~isAFun
        % Command Line Code Approach
        WT = W{1}.hybridFactorVector(t,parameters);
        for i=length(W):-1:1
            % Evaluate the propensity functions at the current state:
            if ~W{i}.isTimeDependent||W{i}.isFactorizable
                props(i+jt) = WT(i)*W{i}.stateDependentFactor(x,parameters);     
            else
                props(i+jt) = W{i}.jointDependentFactor(t,x,parameters);
            end
        end
    else
        % GUI approach
        for i=length(W):-1:1
            % Evaluate the propensity functions at the current state.
            props(i+jt) = W{i}(x,t);     
        end
    end
    
    %% Step 3: Generate uniform random number(s) in [0, 1].
    u1 = rand;
    u2 = rand;

    %% Step 4: Compute the time interval until the next reaction:
    % Note that the sum of the propensities is the inverse of the average
    % waiting time. 

    delta_t = -1 * log(u1) / sum(props);

    %% Step 5: Are there upcoming delayed reactions in [t, t + delta_t]?
    
    scheduledDelayedReactions = sortrows(scheduledDelayedReactions);
    delayedReactionTimes = scheduledDelayedReactions(:, 1);
    rows = find(...
        delayedReactionTimes >= t & delayedReactionTimes <= (t + delta_t));
    rowsize = size(rows);
    if rowsize(1, 1) > 0
        % There IS a delayed reaction in that time interval.
        % Advance t to the time of the next scheduled delayed reaction.

        t = scheduledDelayedReactions(rows, 1);
        rxn = scheduledDelayedReactions(rows, 2);

        % Add the current state to the record of states for all intervening
        % times:
        while (iprint <= Nt) & (t > T_array(iprint))
            X_array(:, iprint) = x;
            iprint = iprint + 1;
        end

        if (t + delta_t) > max(T_array)
            break % Do not simulate beyond the time boundaries.
        end

        % Update the state according to the delayed reaction. First, undo
        % the subtraction of reactants initially performed when the
        % reaction was first scheduled; then apply the stochiometry of the
        % overall reaction.

        x = x - delayedReactionSchedulingS(:, rxn) + S(:, rxn);
    else % rowsize(1, 1) > 0
        %% Step 6: Find the channel of the next reaction:

        threshold = u2 * sum(props); % threshold is in [0, sum(props)]
        rxn = 1;
        while sum(props(1: rxn)) < threshold
            rxn = rxn + 1;
        end
        % At this point, rxn is the number of the chosen reaction.

        %% Step 7: Determine whether the next reaction is delayed:

        rxnIsDelayed = false;
        rxnDelayIndex = 0;
        if hasDelayedReactions
            rxnDelayIndex = find(delayedReactions(:, 1) == rxn, 1);
            rxnIsDelayed = ~isempty(rxnDelayIndex);
        end

        if rxnIsDelayed
            % If the reaction is delayed, postpone the update to
            % t_d = t + tau, where tau is the delay corresponding to the
            % reaction. However, subtract the reactants from the state
            % immediately, so that they cannot be involved in any other
            % reactions in the interim.

            x = x + delayedReactionSchedulingS(:, rxn);
            nextRxn = [t + delayedReactions(rxnDelayIndex, 2), rxn];
            scheduledDelayedReactions = ... 
                [scheduledDelayedReactions; nextRxn];
        else
            % If the reaction is not delayed, perform it; i.e.

            % Advance t to the time of the reaction:

            t = t + delta_t;

            % Add the current state to the record of states for all
            % intervening times:

            while (iprint <= Nt) & (t > T_array(iprint))
                X_array(:, iprint) = x;
                iprint = iprint + 1;
            end

            if t > max(T_array)
                break % Do not simulate beyond the time boundaries.
            end

            % Update the state according to the reaction:

            x = x + S(:, rxn);          
        end % Step 7 ...
    end % rowsize(1, 1) == 0 [No delayed reactions were upcoming]
end % while t < max(T_array) ...