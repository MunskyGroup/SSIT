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
    %       array of size (number of states) x (number of constraints *
    %       number of reactions). The (i,j)-th element equals 1 if
    %       states(:,i) + stoichMatrix(:, k) violates the c-th constraint,
    %       where j = c + (number of constraints - 1)*k.
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
        function obj = FiniteStateSet(states, stoichMatrix)
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
            
            if size(stoichMatrix, 1) ~= size(states, 1)
                error('Stoichiometry matrix and state dimensions mismatch.');
            end
            
            obj.states = states;
            obj.stoichMatrix = stoichMatrix;
            obj.reachableIndices = zeros(size(states,2), size(stoichMatrix, 2));
            key_set = state2key(uint64(states));

            if max([key_set{:}])>1e19
                disp({'WARNING - State index is above machine precision.';'Results may be inaccurate';'Try re-ordering species from low to high expected values'});
            end
            
            obj.state2indMap = containers.Map(key_set, 1:size(states,2));
            if size(obj.states,2)~=obj.state2indMap.Count
                error('HERE')
            end
        end
        
        function obj = expand(obj, fConstraints, bConstraints)
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
            nReactions = size(obj.stoichMatrix, 2);
            obj.numConstraints = nConstraints;
            obj.reachableIndices(obj.reachableIndices<=0) = 0;
            obj.outboundTransitions = zeros(size(obj.states, 2), nReactions*nConstraints, 'uint8');
            
            % We will narrow the search to states reachable from the subset states(:, exploration_range)
            activeNodes = 1:size(obj.states,2);
            
            stop = false;

            if size(obj.states,2)~=obj.state2indMap.Count
                error('Stateset does not match index map.')
            end


            while (~stop)
                n_states_old = size(obj.states, 2);
                
                % search reachable states from current states
                for k = 1:nReactions
                    
                    % We only compute candidates from states that cannot
                    % reach existing states through the current reaction
                    % channel
                    idxToSearch = activeNodes(obj.reachableIndices(activeNodes, k) == 0);

                    candidates = obj.states(:,idxToSearch) + repmat(obj.stoichMatrix(:, k), 1, length(idxToSearch));

                    fVal = fConstraints(candidates);
                    
                    % compute the keys associated with the candidates
                    keySet = state2key(candidates);
                    
                    % check whether the candidate states already exist
                    stateFound = isKey(obj.state2indMap, keySet);
                    stateLocations = cell2mat(values( obj.state2indMap, keySet(stateFound) ));
                    
                    obj.reachableIndices(idxToSearch(stateFound), k) = stateLocations;

                    obj.outboundTransitions(idxToSearch, 1 + (k-1)*nConstraints:k*nConstraints) = uint8(fVal > repmat(bConstraints, 1, size(fVal,2)))';
                    constraintsCheck = (sum(obj.outboundTransitions(idxToSearch, 1 + (k-1)*nConstraints:k*nConstraints), 2)==0)';
                    obj.reachableIndices(idxToSearch(constraintsCheck==0), k) = -1;

                    newStatesCheck = ~stateFound;
                    i_accept_new = find((constraintsCheck == 1) & (newStatesCheck));
                   
                    if (~isempty(i_accept_new))                                                
                         obj.reachableIndices(idxToSearch(i_accept_new), k) = size(obj.states, 2) + (1:length(i_accept_new));
                         
                         for gh = length(i_accept_new):-1:1
                             obj.state2indMap(keySet{i_accept_new(gh)})=(size(obj.states, 2)  + gh);
                         end

                        obj.states = [obj.states candidates(:,i_accept_new)];
                        obj.reachableIndices = [ obj.reachableIndices; 
                                                 zeros(length(i_accept_new), nReactions) ];
                        obj.outboundTransitions = [obj.outboundTransitions; 
                                                   zeros(length(i_accept_new), nConstraints*nReactions)];
                    end

                                        
                end
                if size(obj.states,2)~=obj.state2indMap.Count
                    error('Stateset does not match index map.')
                end
                
                if (size(obj.states, 2) == n_states_old)
                    stop = true;
                else
                    activeNodes = (n_states_old+1: size(obj.states,2));
                end
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

function keys =  state2key( states )
% Hash function to convert a N-dimensional integer vector into a unique
% integer id using recursive Cantor pairing.
%
% References
% ----------
% A. Gupta, J. Mikelson, and M. Khammash, "A finite state projection algorithm for the stationary solution of the chemical master equation," The Journal of Chemical Physics, vol. 147, no. 15, p. 154101, Oct. 2017, doi: 10.1063/1.5006484.
% L. Meri, "Some remarks on the Cantor pairing function,” Le Matematiche, vol. 62, Dec. 2007.
 
d = size(states, 1);
keys = states(1, :);
for i = 2:d
    keys = (keys + states(i, :)).*(keys + states(i, :) + 1)./2 + states(i, :);
end
keys = num2cell(keys);
end


