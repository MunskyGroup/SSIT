classdef FiniteStateSet
    % Object for managing and exploring reachable states of the chemical
    % master equation.
    %
    %
    % Parameters
    % ----------
    %
    %   states: 2-D array
    %       stores the states discovered so far. Each column represent a
    %       state. Thus the number of rows of this array must equal to the
    %       number of species.
    %
    %   stoichMatrix: 2-D array
    %       stores the stoichiometry matrix. Each column correspond to a
    %       stoichiometric vector.
    %
    %   reachableIndices: 2-D array
    %       array of size (number of states) x (number of reactions). The
    %       (i,j) element stores the index of the state reachable by state
    %       i via the j-th reaction.
    %
    %   outboundTransitions: 2-D array
    %       sparse array of size (number of states) x (number of
    %       constraints * number of reactions). The (i,j)-th element
    %       equals 1 if states(:,i) + stoichMatrix(:, k) violates the
    %       c-th constraint, where j = c + (number of constraints - 1)*k.
    %
    %   state2indMap: containers.Map
    %       MATLAB map object for fast state lookup. `state2indMap(x)`
    %       returns the index of state `x` in `states` if exists.
    %
    %   numConstraints: integer
    %       number of shape constraints.
    %
    properties
        states; 
        stoichMatrix;
        reachableIndices; 
        outboundTransitions; 
        state2indMap; 
        numConstraints;                
    end
    
    methods
        function obj = FiniteStateSet(states, stoichMatrix, specialEvents)
            arguments
                states
                stoichMatrix
                specialEvents = []
            end
            % Construct an instance of class FiniteStateSet.
            %
            % Parameters
            % -----------
            %
            % states: 2-D array
            %   array of states, each column represents one state.
            %
            % stoichMatrix: 2-D array
            %   stoichiometry matrix.
            %
            % specialEvents: structure containing options for special
            %   events like geometric bursts or division.
            %
            % Returns
            % -------
            %
            % obj: FiniteStateSet    
            %
            % Examples
            % --------
            % The following code snippet constructs a truncated state space
            % for the two-state bursting gene expression model
            %
            % >>> states = [1 0 0;1 0 1;1 0 2;0 1 0;0 1 1;0 1 2]';
            % >>> stoichMatrix = [-1 1 0;1 -1 0;0 0 1;0 0 -1]';
            % >>> stateSpace = ssit.FiniteStateSet(states, stoichMatrix);         
            
            if ~isempty(stoichMatrix)&&(size(stoichMatrix, 1) ~= size(states, 1))
                error('Stoichiometry matrix and state dimensions mismatch.');
            end
            
            obj.states = states;
            obj.stoichMatrix = stoichMatrix;
            obj.reachableIndices = zeros(size(states,2), size(stoichMatrix, 2)+length(specialEvents), 'int32');
            key_set = num2cell(uint64(states), 1)';

            if max([key_set{:}])>1e19
                disp({'WARNING - State index is above machine precision.';'Results may be inaccurate';'Try re-ordering species from low to high expected values'});
            end
            
            % obj.state2indMap = dictionary(string(key_set), 1:size(states,2));
            obj.state2indMap = dictionary(key_set, (1:size(states,2))');
            % obj.state2indMap = containers.Map(key_set, 1:size(states,2));
            % if size(obj.states,2)~=obj.state2indMap.Count
            if size(obj.states,2)~=obj.state2indMap.numEntries
                error('HERE')
            end
        end
        
        function obj = expand(obj, fConstraints, bConstraints, specialEvents)
            arguments
                obj
                fConstraints
                bConstraints
                specialEvents = [];
            end
            % Expand to all reachable states that satisfy the constraints given by the system
            %       ``f(x) <= b (1)``
            %
            % Parameters
            % ----------
            %
            %   fConstraints: function handle
            %       callable to evaluate the left side of the inequality (1). Must
            %       have syntax
            %           ``Y = f(X)``
            %       where X is an array of states arranged column-wise, and Y
            %       is a matrix, with size(Y,2) == size(X, 2) and size(Y,1) equals number of constraints.
            %
            %   bConstraints: column vector
            %       the right hand side of the inequality (1).
            %
            %   specialEvents: (optional) structure containing parameters and
            %       arguments for special events.
            %
            % Returns
            % -------
            %
            % obj: FiniteStateSet
            %       the expanded version of the calling FiniteStateSet object.
            %
            % Example
            % -------
            %   
            %   
            %   >>> f = @(X) [X(1,:); X(2,:); X(3,:)] ;
            %   >>> b = [10 10 10]';
            %   >>> stateSet = stateSet.expand(f, b);
            %   
            %
            
            nConstraints = numel(bConstraints);
            nRegularRxns = size(obj.stoichMatrix, 2);
            nReactions = nRegularRxns + length(specialEvents);
            obj.numConstraints = nConstraints;
            obj.reachableIndices(obj.reachableIndices<=0) = 0;
            obj.outboundTransitions = sparse(size(obj.states, 2), nReactions*nConstraints);
            nSpecies = size(obj.states,1);
            outboundRows = {};
            outboundCols = {};
            
            % We will narrow the search to states reachable from the subset states(:, exploration_range)
            activeNodes = 1:size(obj.states,2);
            
            stop = false;

            if size(obj.states,2)~=obj.state2indMap.numEntries
                error('Stateset does not match index map.')
            end

            specialModifiers = zeros(1, length(specialEvents));
            for eventIndex = 1:length(specialEvents)
                if ~isfield(specialEvents(eventIndex),'type') || strcmp(specialEvents(eventIndex).type,'decay')
                    specialModifiers(eventIndex) = -1;
                elseif strcmp(specialEvents(eventIndex).type,'production')
                    specialModifiers(eventIndex) = 1;
                else
                    error('Special Event Type not recognized -- must be "production" or "decay"')
                end
            end

            % Detect constraint function mode once (supports two-arg logical output or one-arg numeric)
            useDirectMask = false;
            if ~isempty(obj.states)
                try
                    testMask = fConstraints(obj.states(:,1:1), bConstraints);
                    useDirectMask = islogical(testMask) && size(testMask,1)==nConstraints && size(testMask,2)==1;
                catch
                end
            end

            % Use direct stoichiometry in sequential loop.
            stoich = obj.stoichMatrix;

            while (~stop)
                n_states_old = size(obj.states, 2);

                % Collect candidate transitions to truly new states and add
                % them to the dictionary in one batched operation per iteration.
                allNewCands = zeros(nSpecies, 0, 'like', obj.states);
                allNewRxnIdx = zeros(1, 0);
                allNewK = zeros(1, 0);

                for k = 1:nReactions
                    idxToSearch = activeNodes(obj.reachableIndices(activeNodes, k) == 0);
                    if isempty(idxToSearch)
                        continue
                    end
                    nSearch = numel(idxToSearch);

                    % Build all candidate states for reaction k at once (no inner chunk loop).
                    % For regular reactions, each source state maps to exactly one unique
                    % candidate, so there are no intra-reaction duplicates.
                    if k <= nRegularRxns
                        cands  = obj.states(:, idxToSearch) + stoich(:, k);
                        rxnIdx = idxToSearch;
                    else
                        modifier = specialModifiers(k - nRegularRxns);
                        cands  = zeros(nSpecies, nSearch * nSpecies, 'like', obj.states);
                        rxnIdx = repmat(idxToSearch, 1, nSpecies);
                        for iSp = 1:nSpecies
                            cols = (iSp-1)*nSearch + (1:nSearch);
                            cands(:, cols) = obj.states(:, idxToSearch);
                            cands(iSp, cols) = cands(iSp, cols) + modifier;
                        end
                    end

                    % Evaluate constraints
                    if useDirectMask
                        oMask = fConstraints(cands, bConstraints);
                    else
                        oMask = bsxfun(@gt, fConstraints(cands), bConstraints);
                    end
                    cCheck = reshape(~any(oMask, 1), 1, []);
                    [vCons, vCands] = find(oMask);

                    if ~isempty(vCons)
                        outboundRows{end+1,1} = rxnIdx(vCands(:));
                        outboundRows{end,1} = outboundRows{end,1}(:);
                        outboundCols{end+1,1} = vCons(:) + (k-1)*nConstraints;
                        outboundCols{end,1} = outboundCols{end,1}(:);
                    end

                    % Mark states that lead out of the domain.
                    obj.reachableIndices(rxnIdx(~cCheck), k) = -1;

                    % Key lookup only for candidates that pass constraints.
                    i_feasible = find(cCheck);
                    if ~isempty(i_feasible)
                        feasCands = cands(:, i_feasible);
                        feasRxnIdx = rxnIdx(i_feasible);

                        % Deduplicate feasible states before dictionary lookup.
                        [uniqFeasRows, ~, icFeas] = unique(feasCands', 'rows', 'stable');
                        uniqFeas = uniqFeasRows';
                        kSet = num2cell(uint64(uniqFeas), 1)';
                        sFoundUniq = reshape(isKey(obj.state2indMap, kSet), 1, []);

                        locPerUniq = zeros(1, numel(sFoundUniq), 'int32');
                        if any(sFoundUniq)
                            locPerUniq(sFoundUniq) = int32(obj.state2indMap(kSet(sFoundUniq)));
                        end
                        locPerFeasible = locPerUniq(icFeas);

                        % Set transitions to already existing states now.
                        i_foundFeasible = locPerFeasible > 0;
                        if any(i_foundFeasible)
                            obj.reachableIndices(feasRxnIdx(i_foundFeasible), k) = locPerFeasible(i_foundFeasible);
                        end

                        % Defer insertion of truly new states to a single batched step.
                        i_newFeasible = ~i_foundFeasible;
                        if any(i_newFeasible)
                            allNewCands = [allNewCands, feasCands(:, i_newFeasible)]; %#ok<AGROW>
                            allNewRxnIdx = [allNewRxnIdx, feasRxnIdx(i_newFeasible)]; %#ok<AGROW>
                            allNewK = [allNewK, repmat(k, 1, nnz(i_newFeasible))]; %#ok<AGROW>
                        end
                    end
                end

                % Add all newly discovered states once and map transitions.
                if ~isempty(allNewCands)
                    [uniqRows, ~, icAll] = unique(allNewCands', 'rows', 'stable');
                    uniqCands = uniqRows';
                    nAdd = size(uniqCands, 2);
                    base = size(obj.states, 2);

                    obj.states = [obj.states, uniqCands];
                    obj.reachableIndices = [obj.reachableIndices; zeros(nAdd, nReactions, 'int32')];

                    newIdxVals = int32(base + (1:nAdd));
                    newKeys = num2cell(uint64(uniqCands), 1)';
                    obj.state2indMap(newKeys) = newIdxVals(:);

                    nStatesNow = size(obj.reachableIndices, 1);
                    linIdx = allNewRxnIdx(:) + (allNewK(:)-1) * nStatesNow;
                    obj.reachableIndices(linIdx) = newIdxVals(icAll(:));
                end

                if size(obj.states,2)~=obj.state2indMap.numEntries
                    error('Stateset does not match index map.')
                end

                if (size(obj.states, 2) == n_states_old)
                    stop = true;
                else
                    activeNodes = (n_states_old+1: size(obj.states,2));
                end
            end

            if isempty(outboundRows)
                obj.outboundTransitions = sparse(size(obj.states, 2), nReactions*nConstraints);
            else
                rowIdx = double(vertcat(outboundRows{:}));
                colIdx = double(vertcat(outboundCols{:}));
                obj.outboundTransitions = spones(sparse(rowIdx, colIdx, 1, size(obj.states, 2), nReactions*nConstraints));
            end
        end
                
        function n = getNumStates(obj)
            n = size(obj.states, 2);
        end
        
        function d = getNumSpecies(obj)
            d = size(obj.states, 1);
        end
    end
end