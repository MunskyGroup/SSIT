function [X_array] = runSingleSsaModifiedTwo(...
    x0, S, W, T_array, isTimeVarying, ...
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
% time; the second, the reaction number. Initially, the array will contain
% a million rows; we will grow it over time. Note that the correct default
% values are inf for time and zero for reaction number (since we use
% one-based indexing).

if hasDelayedReactions
    sdrRowCount = 1e6;
    % This is the first row that is valid for scheduling new reactions.
    sdrInsertRow = 1;
else
    sdrRowCount = 0;
    sdrInsertRow = 0;
end

scheduledDelayedReactions = ones(sdrRowCount, 2);
scheduledDelayedReactions(:, 1) = inf * scheduledDelayedReactions(:, 1);
scheduledDelayedReactions(:, 2) = 0 * scheduledDelayedReactions(:, 2);

sdrLastAccessedRow = sdrRowCount + 1; % Necessarily invalid

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

    %% Step 5: Find the channel of the next reaction:

    threshold = u2 * sum(props); % threshold is in [0, sum(props)]
        rxn = 1;
        while sum(props(1: rxn)) < threshold
            rxn = rxn + 1;
        end
        % At this point, rxn is the number of the chosen reaction.

    %% Step 6: Determine whether the next reaction is delayed:

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
        % nextRxn = [t + delayedReactions(rxnDelayIndex, 2), rxn];
        scheduledDelayedReactions(sdrInsertRow, 1) = ...
            t + delayedReactions(rxnDelayIndex, 2);
        scheduledDelayedReactions(sdrInsertRow, 2) = rxn;
        sdrInsertRow = sdrInsertRow + 1;
        if sdrInsertRow > sdrRowCount
            if sdrLastAccessedRow > sdrRowCount
                error('Ran out of rows for queueing delayed reactions!')
            else
                % Reinitialize the array and start over.
                scheduledDelayedReactions = ones(sdrRowCount, 2);
                scheduledDelayedReactions(:, 1) = ...
                    inf * scheduledDelayedReactions(:, 1);
                scheduledDelayedReactions(:, 2) = ...
                    0 * scheduledDelayedReactions(:, 2);
                scheduledDelayedReactions(1, 1) = ...
                    t + delayedReactions(rxnDelayIndex, 2);
                scheduledDelayedReactions(1, 2) = rxn;
                sdrInsertRow = 2;
            end
        end
        % scheduledDelayedReactions = ... 
        %     [scheduledDelayedReactions; nextRxn];
    else
        % If the reaction is not delayed, determine whether there is a
        % delayed reaction previously scheduled to occur within the
        % relevant time interval. If there is, the first such delayed
        % reaction is performed; otherwise, this newly selected
        % reaction is performed.

        delayedReactionTimes = scheduledDelayedReactions(...
            sdrLastAccessedRow + 1 : sdrInsertRow, ...
            1);
        
        % This call will find the first row with a time in 
        % [t, t + delta_t], if one exists. As the matrix has been 
        % sorted by time, this will be the earliest such time in the 
        % interval.
    
        row = find(...
            delayedReactionTimes >= t & ...
            delayedReactionTimes <= (t + delta_t), ...
            1);
        rowsize = size(row);
        if rowsize(1, 1) > 0
            % There IS a delayed reaction in that time interval.
            % Advance t to the time of the next scheduled 
            % delayed reaction.
    
            % Shift the row to reflect the earlier 
            % subsetting of the array.

            row = row + sdrLastAccessedRow;
            sdrLastAccessedRow = row;
            t = scheduledDelayedReactions(row, 1);
            rxn = scheduledDelayedReactions(row, 2);
    
            % Having extracted the next scheduled delayed reaction,
            % remove it and all earlier rows from the record.
    
            % Add the current state to the record of states for all
            % intervening times:

            while (iprint <= Nt) & (t > T_array(iprint))
                X_array(:, iprint) = x;
                iprint = iprint + 1;
            end
    
            if t > max(T_array)
                break % Do not simulate beyond the time boundaries.
            end
    
            % Update the state according to the delayed reaction.
            % First, undo the subtraction of reactants initially 
            % performed when the reaction was first scheduled; then 
            % apply the stochiometry of the overall reaction.
    
            x = x - delayedReactionSchedulingS(:, rxn) + S(:, rxn);
        else % rowsize(1, 1) > 0
            % No delayed reactions are upcoming, so perform the
            % newly selected reaction:

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
        end % rowsize(1, 1) == 0 [No delayed reactions were upcoming]            
    end % Step 6 ...  
end % while t < max(T_array) ...