classdef FspMatrixTerm
    % Data structure to represent a term of the truncated CME matrix in the
    % affine decomposition (see documentation of
    % :mat:class:`~+ssit.@FspMatrix.FspMatrix` for details).
    %
    % Parameters
    % ----------
    %
    % isTimeDependent: true/false
    %   Does the matrix depend on time.
    %
    % isFactorizable: true/false
    %   If the matrix is time-dependent, does it factorize into the form
    %   \(c(t)\cdot B\) where \(c(t)\) is a scalar function of time and
    %   \(B\) is a time-invariant matrix.
    %
    % propensity: :mat:class:`~+ssit.@Propensity.Propensity`
    %   The propensity used to generate entries of this matrix.
    %
    % matrix:
    %   Generated sparse matrix.
    %
    % numConstraints: integer
    %    Number of constraints/sink states.
    properties
        isTimeDependent logical
        isFactorizable logical
        propensity ssit.Propensity
        matrix
        numConstraints {mustBeInteger(numConstraints)}
    end

    methods
        function obj = FspMatrixTerm(propensity, stateSet, numConstraints)
            % Construct a `FspMatrixTerm` object from a propensity
            % function, a FSP-truncated state space, and number of
            % constraints..
            %
            % Parameters
            % ----------
            %
            % propensity: :mat:class:`~+ssit.@Propensity.Propensity`
            %   Propensity function of the reaction.
            %
            % stateSet: :mat:class:`~+ssit.@FiniteStateSet.FiniteStateSet`
            %   State space from which to construct the matrix entries.
            %
            % numConstraints: integer
            %   Number of FSP shape constraints or sink states.
            obj.propensity = propensity;
            obj.isTimeDependent = propensity.isTimeDependent;
            obj.isFactorizable = propensity.isFactorizable;
            obj.numConstraints = numConstraints;

            if ((~obj.isTimeDependent) || (obj.isFactorizable))
                obj.matrix = ssit.FspMatrixTerm.GenerateMatrixTermTimeInvariant(propensity, stateSet, numConstraints); % The matrix data is generated only once
            else
                obj.matrix = stateSet; % The data has to be generated with every matrix-vector multiplication
            end
        end

        function w = multiply(obj, t, v)
            % Compute the action of the term matrix at time `t` on a vector
            % `v`.
            if (~obj.isTimeDependent)
                w = obj.matrix*v;
            elseif (obj.isFactorizable)
                w = obj.propensity.timeDependentFactor(t)*(obj.matrix*v);
            else
%                 try
                    A = ssit.FspMatrixTerm.generateTimeVaryingMatrixTerm(t, obj.propensity, obj.matrix, obj.numConstraints);
%                 catch
%                     1+1
%                 end
                w = A*v;
            end
        end
    end

    methods (Static)
        function A_fsp = GenerateMatrixTermTimeInvariant(propensity, state_set, numConstraints)

            numConstraints = state_set.numConstraints;
            n_states = size(state_set.states, 2);
            reachableIndices = state_set.reachableIndices(:, propensity.reactionIndex);

            prop_val = propensity.evaluateStateFactor(state_set.states)...
                +0*zeros(1,size(state_set.states,2));
            if (size(prop_val, 2) > 1)
                prop_val = prop_val';
            end

            % Filling in values for MATLAB's sparse matrix format (which is coo)
            i = zeros(n_states*2, 1);
            j = i;
            aij = zeros(length(i), 1);

            % Coordinates and values for diagonal entries
            j(1:n_states) = 1:n_states;
            i(1:n_states) = 1:n_states;
            aij(1:n_states) = -1.0*prop_val;

            % Coordinates and values for normal states
            offset = n_states;
            j((1:n_states) + offset) = 1:n_states;
            i((1:n_states) + offset) = reachableIndices(1:n_states);
            aij((1:n_states) + offset) = prop_val;

            % Coordinates and values for sink states
            ireaction = propensity.reactionIndex;
            n_constrs_failed = sum(state_set.outboundTransitions(:, 1 + numConstraints*(ireaction-1):numConstraints*ireaction), 2);
            isinks = zeros(nnz(state_set.outboundTransitions(:, 1 + numConstraints*(ireaction-1):numConstraints*ireaction)), 1);
            jsinks = isinks;
            aijsinks = zeros(length(jsinks), 1);
            k = 1;
            for c = 1:state_set.numConstraints
                ninsert = nnz(state_set.outboundTransitions(:, c + numConstraints*(ireaction-1)));
                isinks(k:k+ninsert-1) = n_states + c;
                jsinks(k:k+ninsert-1) = find(state_set.outboundTransitions(:,c + numConstraints*(ireaction-1)));
                try
                    aijsinks(k:k+ninsert-1) = prop_val(jsinks(k:k+ninsert-1))./n_constrs_failed(jsinks(k:k+ninsert-1));
                catch
                    aaa=1;
                end
                k = k + ninsert;

            end

            % eliminate out-of-bound coordinates
            indices_keep = find( (i > 0) & ( j > 0) );
            i = i(indices_keep);
            j = j(indices_keep);
            aij = aij(indices_keep);

            A_fsp = sparse([i;isinks], [j;jsinks], [aij;aijsinks], n_states + numConstraints, n_states + numConstraints);
        end

        function A_fsp = generateTimeVaryingMatrixTerm(t, propensity, state_set, numConstraints)

            n_states = size(state_set.states, 2);
            reachableIndices = state_set.reachableIndices(:, propensity.reactionIndex);
            reachableIndices(reachableIndices < 0) = 0;

            prop_val = propensity.evaluate(t, state_set.states);
            if (size(prop_val, 2) > 1)
                prop_val = prop_val';
            end

            % Filling in values for MATLAB's sparse matrix format (which is coo)
            i = zeros(n_states*2, 1);
            j = i;
            aij = zeros(length(i), 1);

            % Coordinates and values for diagonal entries
            j(1:n_states) = 1:n_states;
            i(1:n_states) = 1:n_states;
            aij(1:n_states) = -1.0*prop_val;

            offset = n_states;
            j((1:n_states) + offset) = 1:n_states;
            i((1:n_states) + offset) = reachableIndices(1:n_states);
            aij((1:n_states) + offset) = prop_val(:);

            % Coordinates and values for sink states
            ireaction = propensity.reactionIndex;
            n_constrs_failed = sum(state_set.outboundTransitions(:, 1 + numConstraints*(ireaction-1):numConstraints*ireaction), 2);
            isinks = zeros(nnz(state_set.outboundTransitions(:, 1 + state_set.numConstraints*(ireaction-1):state_set.numConstraints*ireaction-1)), 1);
            jsinks = isinks;
            aijsinks = zeros(length(jsinks), 1);
            k = 1;
            for c = 1:state_set.numConstraints
                ninsert = nnz(state_set.outboundTransitions(:, c + state_set.numConstraints*(ireaction-1)));
                isinks(k:k+ninsert-1) = n_states + c;
                jsinks(k:k+ninsert-1) = find(state_set.outboundTransitions(:,c + state_set.numConstraints*(ireaction-1)));
                aijsinks(k:k+ninsert-1) = prop_val(jsinks(k:k+ninsert-1))./n_constrs_failed(jsinks(k:k+ninsert-1));

                k = k + ninsert;
            end

            % eliminate out-of-bound coordinates
            indices_keep = find( (i > 0) & ( j > 0) );
            i = i(indices_keep);
            j = j(indices_keep);
            aij = aij(indices_keep);

%             try
                A_fsp = sparse([i;isinks], [j;jsinks], [aij;aijsinks], n_states + numConstraints, n_states + numConstraints);
%             catch
%                 1+1
%             end
            end
    end
end
