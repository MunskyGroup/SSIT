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
    %
    % modRedTransformMatrices: structure
    %   The transformation matrices needed to create a projected (reduced)
    %   model:
    %      A' = PHIinv*A*PHI.*phiScale
    %   see [https://ieeexplore.ieee.org/abstract/document/6425828] for
    %   model reduction details 

    properties
        isTimeDependent logical
        isFactorizable logical
        propensity ssit.Propensity
        matrix
        numConstraints {mustBeInteger(numConstraints)}
        modRedTransformMatrices
    end

    methods
        function obj = FspMatrixTerm(propensity, parameters, stateSet, numConstraints, varNames, modRedTransformMatrices)
            arguments
                propensity
                parameters
                stateSet
                numConstraints
                varNames = []
                modRedTransformMatrices = []
            end
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
            %
            % varNames: cell of strings
            %   list of names of all species in the model.
            %
            % modRedTransformMatrices: structure
            %   The transformation matrices needed to create a projected (reduced)
            %   model:
            %      A' = PHIinv*A*PHI.*phiScale
            %   see [https://ieeexplore.ieee.org/abstract/document/6425828] for
            %   model reduction details
            obj.propensity = propensity;
            obj.isTimeDependent = propensity.isTimeDependent;
            obj.isFactorizable = propensity.isFactorizable;
            obj.numConstraints = numConstraints;
            obj.modRedTransformMatrices = modRedTransformMatrices;

            if ((~obj.isTimeDependent) || (obj.isFactorizable))
                % Generate generator matrix for this term
                obj.matrix = ssit.FspMatrixTerm.GenerateMatrixTermTimeInvariant(propensity, parameters, stateSet, numConstraints, varNames); % The matrix data is generated only once
                
                % Perform model reduction if requested
                if ~isempty(modRedTransformMatrices)
                    obj.matrix = modRedTransformMatrices.phi_inv*...
                        obj.matrix(1:end-numConstraints,1:end-numConstraints)*...
                        modRedTransformMatrices.phi;
                    if ~isempty(obj.modRedTransformMatrices.phiScale)
                        obj.matrix = obj.matrix - diag(diag(obj.matrix));
                        obj.matrix = obj.matrix.*obj.modRedTransformMatrices.phiScale;
                        obj.matrix = obj.matrix - diag(sum(obj.matrix));
                    end
                end
            else
                obj.matrix = stateSet; % The data has to be generated with every matrix-vector multiplication
            end
        end

        function w = multiply(obj, t, v, parameters)
            % Compute the action of the term matrix at time `t` on a vector
            % `v`.
            if (~obj.isTimeDependent)&&isempty(obj.propensity.timeDependentFactor)
                w = obj.matrix*v;
            elseif (~obj.isTimeDependent)&&~isempty(obj.propensity.timeDependentFactor)
                w = obj.propensity.timeDependentFactor(t,parameters)*(obj.matrix*v);
            elseif (obj.isFactorizable)
                w = obj.propensity.timeDependentFactor(t,parameters)*(obj.matrix*v);
            else
                A = ssit.FspMatrixTerm.generateTimeVaryingMatrixTerm(t, obj.propensity, obj.matrix, parameters, obj.numConstraints);
                if ~isempty(obj.modRedTransformMatrices)
                    A = obj.modRedTransformMatrices.phi_inv*...
                        A(1:end-obj.numConstraints,1:end-obj.numConstraints)*...
                        obj.modRedTransformMatrices.phi;
                    if ~isempty(obj.modRedTransformMatrices.phiScale)
                        A = A - diag(diag(A));
                        A = A.*obj.modRedTransformMatrices.phiScale;
                        A = A - diag(sum(A));
                    end
                end
                w = A*v;
            end
        end

        function w = multiplyHybrid(obj, t, v, parameters, vODEs)
            % Compute the action of the term matrix at time `t` on a vector
            % `v`. This version is for hybrid models, where 'v' corresponds
            % to the CME (i.e., the vector P(t)) and'vODEs' corresponds to
            % the upstream ODEs.
            
            if (obj.isFactorizable)
                w1 = obj.propensity.hybridFactor(t,parameters,vODEs)*(obj.matrix*v);
            else
                A = ssit.FspMatrixTerm.generateHybridMatrixTerm(t, obj.propensity, obj.matrix, parameters, obj.numConstraints, vODEs');
                if ~isempty(obj.modRedTransformMatrices)
                    A = obj.modRedTransformMatrices.phi_inv*...
                        A(1:end-obj.numConstraints,1:end-obj.numConstraints)*...
                        obj.modRedTransformMatrices.phi;
                    if ~isempty(obj.modRedTransformMatrices.phiScale)
                        A = A - diag(diag(A));
                        A = A.*obj.modRedTransformMatrices.phiScale;
                        A = A - diag(sum(A));
                    end
                end
                w1 = A*v;
            end
            if max(abs(obj.propensity.ODEstoichVector))~=0
                w2 = obj.propensity.ODEstoichVector*obj.propensity.hybridFactor(t,parameters,vODEs);
            else
                w2 = 0*vODEs';
            end
            w = [w1;w2];

        end
    end

methods (Static)

    function A_fsp = GenerateMatrixTermTimeInvariant(propensity, parameters, state_set, numConstraints, varNames)
        arguments
            propensity
            parameters
            state_set
            numConstraints =[];
            varNames = [];
        end
        numConstraints = state_set.numConstraints;
        n_states = size(state_set.states, 2);
        reachableIndices = state_set.reachableIndices(:, propensity.reactionIndex);

%         if propensity.isTimeDependent&&~isempty(propensity.timeDependentFactor)
%             % Constants are sometimes lumped into a constant time dependent
%             % term for later convenience.
%             prop_val = propensity.timeDependentFactor(0).*...
%                 propensity.evaluateStateFactor(state_set.states,varNames)...
%                 +0*zeros(1,size(state_set.states,2));
%         else
            prop_val = propensity.evaluateStateFactor(state_set.states,parameters,varNames)...
                +0*zeros(1,size(state_set.states,2));
%         end

        if min(prop_val)<0
            warning('PropFun:negProp','Negative propensity function detected. Results may be incorrect.');
            warning('OFF','PropFun:negProp')
            prop_val(prop_val<0) = 0;
        end

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
                aijsinks(k:k+ninsert-1) = prop_val(jsinks(k:k+ninsert-1))./n_constrs_failed(jsinks(k:k+ninsert-1));
            k = k + ninsert;
        end

        % eliminate out-of-bound coordinates
        indices_keep = find( (i > 0) & ( j > 0) );
        i = i(indices_keep);
        j = j(indices_keep);
        aij = aij(indices_keep);

        A_fsp = sparse([i;isinks], [j;jsinks], [aij;aijsinks], n_states + numConstraints, n_states + numConstraints);
    end

    function A_fsp = generateTimeVaryingMatrixTerm(t, propensity, state_set, parameters, numConstraints, modRedTransformMatrices)
        arguments
            t
            propensity
            state_set
            parameters
            numConstraints
            modRedTransformMatrices =[];
        end

        n_states = size(state_set.states, 2);
        reachableIndices = state_set.reachableIndices(:, propensity.reactionIndex);
        reachableIndices(reachableIndices < 0) = 0;

        prop_val = propensity.evaluate(t, state_set.states, parameters);
        if (size(prop_val, 2) > 1)
            prop_val = prop_val';
        end
        if size(prop_val,1)==1&&size(state_set.states,2)>1
            prop_val = prop_val*ones(size(state_set.states,2),1);
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

        A_fsp = sparse([i;isinks], [j;jsinks], [aij;aijsinks], n_states + numConstraints, n_states + numConstraints);
        if ~isempty(modRedTransformMatrices)
            A_fsp = modRedTransformMatrices.phi_inv*...
                A_fsp(1:end-numConstraints,1:end-numConstraints)*...
                modRedTransformMatrices.phi;
            if ~isempty(modRedTransformMatrices.phiScale)
                A_fsp = A_fsp - diag(diag(A_fsp));
                A_fsp = A_fsp.*modRedTransformMatrices.phiScale;
                A_fsp = A_fsp - diag(sum(A_fsp));
            end
        end


    end

    function A_fsp = generateHybridMatrixTerm(t, hybridPropensity, state_set, parameters, numConstraints, upStreamSpecies, modRedTransformMatrices)
        arguments
            t
            hybridPropensity
            state_set
            parameters
            numConstraints
            upStreamSpecies
            modRedTransformMatrices =[];
        end

        n_states = size(state_set.states, 2);
        reachableIndices = state_set.reachableIndices(:, hybridPropensity.reactionIndex);
        reachableIndices(reachableIndices < 0) = 0;

        prop_val = hybridPropensity.evaluate(t, state_set.states, parameters, upStreamSpecies);
        if (size(prop_val, 2) > 1)
            prop_val = prop_val';
        end
        if size(prop_val,1)==1&&size(state_set.states,2)>1
            prop_val = prop_val*ones(size(state_set.states,2),1);
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
        ireaction = hybridPropensity.reactionIndex;
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

        A_fsp = sparse([i;isinks], [j;jsinks], [aij;aijsinks], n_states + numConstraints, n_states + numConstraints);
        if ~isempty(modRedTransformMatrices)
            A_fsp = modRedTransformMatrices.phi_inv*...
                A_fsp(1:end-numConstraints,1:end-numConstraints)*...
                modRedTransformMatrices.phi;
            if ~isempty(modRedTransformMatrices.phiScale)
                A_fsp = A_fsp - diag(diag(A_fsp));
                A_fsp = A_fsp.*modRedTransformMatrices.phiScale;
                A_fsp = A_fsp - diag(sum(A_fsp));
            end
        end

    end

end
end
