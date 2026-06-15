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
            obj.state2indMap = dictionary(key_set, 1:size(states,2));
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

            regularChunkSize = max(1, floor(1e6 / max(nSpecies, 1)));
            specialChunkSize = max(1, floor(1e6 / max(nSpecies * nSpecies, 1)));
            useDirectMask = [];


            while (~stop)
                n_states_old = size(obj.states, 2);
                
                % search reachable states from current states
                for k = 1:nReactions
                    % We only compute candidates from states that cannot
                    % reach existing states through the current reaction
                    % channel
                    idxToSearch = activeNodes(obj.reachableIndices(activeNodes, k) == 0);
                    if isempty(idxToSearch)
                        continue
                    end

                    if k<=nRegularRxns
                        chunkSize = regularChunkSize;
                    else
                        chunkSize = specialChunkSize;
                    end

                    pendingNewStates = zeros(nSpecies, 0, 'like', obj.states);
                    pendingCount = 0;

                    for idxStart = 1:chunkSize:length(idxToSearch)
                        idxEnd = min(idxStart + chunkSize - 1, length(idxToSearch));
                        currentIdx = idxToSearch(idxStart:idxEnd);
                        nCurrent = length(currentIdx);

                        if k<=nRegularRxns
                            candidates = bsxfun(@plus, obj.states(:,currentIdx), obj.stoichMatrix(:, k));
                            reactionIdx = currentIdx;
                        else
                            modifier = specialModifiers(k-nRegularRxns);
                            candidates = zeros(nSpecies, nCurrent*nSpecies, 'like', obj.states);
                            reactionIdx = repmat(currentIdx, 1, nSpecies);
                            currentStates = obj.states(:,currentIdx);
                            for iSp = 1:nSpecies
                                colRange = (iSp-1)*nCurrent + (1:nCurrent);
                                candidates(:,colRange) = currentStates;
                                candidates(iSp,colRange) = candidates(iSp,colRange) + modifier;
                            end
                        end

                        if isempty(useDirectMask)
                            testCols = min(1, size(candidates,2));
                            try
                                testMask = fConstraints(candidates(:,1:testCols), bConstraints);
                                useDirectMask = islogical(testMask) && isequal(size(testMask), [nConstraints, testCols]);
                            catch
                                useDirectMask = false;
                            end
                        end

                        if useDirectMask
                            outboundMask = fConstraints(candidates, bConstraints);
                        else
                            fVal = fConstraints(candidates);
                            outboundMask = bsxfun(@gt, fVal, bConstraints);
                        end
                        [violatedConstraints, violatedCandidates] = find(outboundMask);
                        if ~isempty(violatedConstraints)
                            outboundRows{end+1,1} = reactionIdx(violatedCandidates(:));
                            outboundRows{end,1} = outboundRows{end,1}(:);
                            outboundCols{end+1,1} = violatedConstraints(:) + (k-1)*nConstraints;
                            outboundCols{end,1} = outboundCols{end,1}(:);
                        end
                        constraintsCheck = reshape(~any(outboundMask, 1), 1, []);

                        % compute the keys associated with the candidates
                        keySet = num2cell(uint64(candidates), 1)';

                        % check whether the candidate states already exist
                        stateFound = isKey(obj.state2indMap, keySet);
                        stateFound = reshape(stateFound, 1, []);
                        if any(stateFound)
                            stateLocations = obj.state2indMap(keySet(stateFound));
                            obj.reachableIndices(reactionIdx(stateFound), k) = stateLocations;
                        end
                        
                        obj.reachableIndices(reactionIdx(~constraintsCheck), k) = -1;

                        newStatesCheck = ~stateFound;
                        i_accept_new = find(constraintsCheck & newStatesCheck);

                        % Remove duplicates if any.
                        if ~isempty(i_accept_new)
                            [~,ia] = unique(candidates(:,i_accept_new)','stable','rows');
                            i_accept_new = i_accept_new(ia);
                        end

                        if (~isempty(i_accept_new))
                             newStateStart = size(obj.states, 2) + pendingCount;
                             obj.reachableIndices(reactionIdx(i_accept_new), k) = newStateStart + (1:length(i_accept_new));
                             for gh = length(i_accept_new):-1:1
                                 obj.state2indMap(keySet(i_accept_new(gh)))=(newStateStart  + gh);
                             end

                            pendingBlock = candidates(:,i_accept_new);
                            pendingNewStates = [pendingNewStates pendingBlock];
                            pendingCount = pendingCount + size(pendingBlock, 2);
                        end
                    end

                    if pendingCount > 0
                        obj.states = [obj.states pendingNewStates];
                        obj.reachableIndices = [ obj.reachableIndices;
                                                 zeros(pendingCount, nReactions, 'int32') ];
                    end
                end
                % if size(obj.states,2)~=obj.state2indMap.Count
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