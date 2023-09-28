classdef FspMatrix
    % The FSP-truncated transition rate matrix of the CME. This data
    % structure is based on decomposing the CME matrix in the form
    %
    % .. math::
    %   A(t, \theta) = A_1(t, \theta) + \ldots + A_R(t, \theta)
    %
    % where \\(t\\) is the time variable, \\(\theta\\) the vector of model parameters, and
    % \\(R\\) is the numer of reactions.
    %
    % Parameters
    % ----------
    %
    % terms: cell array of :mat:class:`~+ssit.@FspMatrixTerm.FspMatrixTerm` instances
    %   FSP matrix terms. Each term correspond to a reaction.
    %
    % See also
    % --------
    %
    % :mat:class:`~+ssit.@FspMatrixTerm.FspMatrixTerm`
    properties
        terms;
    end

    methods
        function obj = FspMatrix(propensities, parameters, stateSet, numConstraints, varNames, modRedTransformMatrices)
            % Construct an instance of FspMatrix.
            %
            % Parameters
            % ----------
            %
            %   propensities: cell of :mat:class:`~+ssit.@Propensity.Propensity` objects
            %       cell of propensities that define the stochastic reaction network
            %       model.
            %
            %   stateSet: an instance of class :mat:class:`~+ssit.@FiniteStateSet.FiniteStateSet`.
            %       the states based on which this transition rate matrix is
            %       built.
            %
            %   numConstraints: integer
            %       number of FSP constraints.
            %
            %   varNames: cell of strings
            %       names of the species considered in the model
            %
            %   modRedTransformMatrices: structure
            %       optional structure containing the projection
            %       transformation matrices PHI and PHIinv to be use for
            %       model reduction.
            %
            % Returns
            % -------
            %
            %   obj: an instance of this class.
            %
            arguments
                propensities
                parameters
                stateSet
                numConstraints
                varNames = []
                modRedTransformMatrices =[];
            end
            obj = obj.regenerate(propensities, parameters, stateSet, numConstraints, varNames, modRedTransformMatrices);
        end

        function obj = regenerate(obj, propensities, stateSet, parameters, numConstraints, varNames, modRedTransformMatrices)
            % Regenerate the state space for the current contraints.
            % Parameters
            % ----------
            %
            %   propensities: cell of :mat:class:`~+ssit.@Propensity.Propensity` objects
            %       cell of propensities that define the stochastic reaction network
            %       model.
            %
            %   stateSet: an instance of class :mat:class:`~+ssit.@FiniteStateSet.FiniteStateSet`.
            %       the states based on which this transition rate matrix is
            %       built.
            %
            %   numConstraints: integer
            %       number of FSP constraints.
            %
            %   varNames: cell of strings
            %       names of the species considered in the model
            %
            %   modRedTransformMatrices: structure
            %       optional structure containing the projection
            %       transformation matrices PHI and PHIinv to be use for
            %       model reduction.
            %
            % Returns
            % -------
            %
            %   obj: an instance of this class.
            arguments
                obj
                propensities
                stateSet
                parameters
                numConstraints
                varNames = []
                modRedTransformMatrices = []
            end
            n_reactions = length(propensities);
            obj.terms = cell(n_reactions, 1);
            for i = 1:n_reactions
                obj.terms{i} = ssit.FspMatrixTerm(propensities{i}, stateSet, parameters, numConstraints, varNames, modRedTransformMatrices);
            end
        end

        function w = multiply(obj,t,v,parameters)
            % Compute the action of the time-dependent FSP matrix.
            %
            % Parameters
            % ----------
            %
            %   t: double
            %       time to evaluate the matrix-vector product.
            %
            %   v: column vector
            %       the input vector.
            %
            % Returns
            % -------
            %
            %   w: column vector.
            %       the output vector.
            wt = obj.terms{1}.propensity.hybridFactorVector(t,parameters);
            w = 0*v;
            % w = obj.terms{1}.multiply(t,v,parameters);
            for i = 1:length(obj.terms)
                if obj.terms{i}.isTimeDependent
                    w = w + obj.terms{i}.multiply(t,v,parameters);
                else
                    w = w + wt(i)*(obj.terms{i}.matrix*v);
                end
            end
        end

        function w = hybridRHS(obj,t,v,parameters,upstreamODEs)
            % Compute the action of the time-dependent FSP matrix.
            %
            % Parameters
            % ----------
            %
            %   t: double
            %       time to evaluate the matrix-vector product.
            %
            %   v: column vector
            %       the input vector.
            %
            %   upstreamODEs: cell of strings
            %       names of the upstream species that will be treated as
            %       ODES. 
            %
            % Returns
            % -------
            %
            %   w: column vector.
            %       the output vector.

            wt = obj.terms{1}.propensity.hybridFactorVector(t,parameters,v(end-length(upstreamODEs)+1:end)');
            vJ1 = v(1:length(v)-length(upstreamODEs));
            w = wt(1)*[obj.terms{1}.matrix*vJ1;...
                obj.terms{1}.propensity.ODEstoichVector];
            for i = 2:length(obj.terms)
                w = w + wt(i)*[obj.terms{i}.matrix*vJ1;...
                    obj.terms{i}.propensity.ODEstoichVector];
            end

        end

        function A = createSingleMatrix(obj, t, parameters, modRedTransformMatrices)
            % Generate a single MATLAB sparse matrix from the current time
            % `t`. This matrix's left multiplication on a vector `v` of
            % appropriate size will return the same output as when calling
            % `multiply(obj,t, v)`.
            % if modRedTransformMatrices is provided, then the model
            % reduction transofrmation will be carried out to return the
            % reduce system.
            arguments
                obj
                t
                parameters
                modRedTransformMatrices =[];
            end

            if obj.terms{1}.isFactorizable
                if (obj.terms{1}.isTimeDependent)
                    wt = obj.terms{1}.propensity.hybridFactorVector(t,parameters);
                    A = wt(1)*obj.terms{1}.matrix;
                else
                    wt = obj.terms{1}.propensity.hybridFactorVector(0,parameters);
                    A = wt(1)*obj.terms{1}.matrix;
                end
            else
                A = ssit.FspMatrixTerm.generateTimeVaryingMatrixTerm(t, obj.terms{1}.propensity, obj.terms{1}.matrix, parameters, obj.terms{1}.numConstraints, modRedTransformMatrices);
            end

            for i = 2:length(obj.terms)
                if obj.terms{i}.isFactorizable
                    % if obj.terms{i}.isTimeDependent
                        A = A + wt(i)*obj.terms{i}.matrix;
                    % else
                        % A = A + obj.terms{i}.propensity.timeDependentFactor(0,parameters)*obj.terms{i}.matrix;
%                         A = A + obj.terms{i}.matrix;
                    % end
                else
                    A = A + ssit.FspMatrixTerm.generateTimeVaryingMatrixTerm(t, obj.terms{i}.propensity, obj.terms{i}.matrix, parameters, obj.terms{i}.numConstraints, modRedTransformMatrices);
                end
            end
            
            ANaN = isnan(A);
            if sum(ANaN(:))>0
                warning('NaNs detected and set to zero. Errors may be present.')
                A(ANaN) = 0;
            end
            
        end
    
        function A = createJacHybridMatrix(obj, t, v, parameters, nUpstream, partialJacobian)
            % Generate a single MATLAB sparse matrix from the current time
            % `t`. This matrix's left multiplication on a vector `v` of
            % appropriate size will return the same output as when calling
            % `multiply(obj,t, v)`.
            arguments
                obj
                t
                v
                parameters
                nUpstream
                partialJacobian = false;
            end

            if partialJacobian
                A = hybridMatrix(obj,t,parameters,v);
                return
            else
                v1 = v(1:end-nUpstream);
                v2 = v(end-nUpstream+1:end);
                Sode = zeros(nUpstream,length(obj.terms));
                dwdv = obj.terms{1}.propensity.DhybridFactorDodesVec(t,parameters,v2');
                gB = sparse(zeros(length(v1),length(v2)));
                for i = 1:length(obj.terms)
                    Sode(:,i) = obj.terms{i}.propensity.ODEstoichVector;
                    if ~max(abs(obj.terms{i}.propensity.stoichVector)>0)
                        gB = gB + sparse(obj.terms{i}.matrix*v1*dwdv(i,:));
                    end
                end
                gD = Sode*dwdv;

                gA = hybridMatrix(obj,t,parameters,v2);
                gC = zeros(length(v2),length(v1));
                % gB = zeros(length(v1),length(v2));
                % for i = 1:length(v2)
                    % delt = max(1e-6,abs(v2(i))/1000);
                    % v2p = v2+delt;
                    % gB(:,i) = (hybridMatrix(obj,t,parameters,v2p)-gA)*v1/delt;

                    % w = obj.terms{1}.propensity.hybridFactorVector(t,parameters,v2');
                    % w2 = obj.terms{1}.propensity.hybridFactorVector(t,parameters,v2p');
                    % wp = (w2-w)/delt;
                    
                    % for k = 1:length(obj.terms)
                    %     % gD(:,i) =  gD(:,i) + obj.terms{k}.propensity.ODEstoichVector*...
                    %     %     (obj.terms{k}.propensity.hybridFactor(t,v2p)-...
                    %     %     obj.terms{k}.propensity.hybridFactor(t,v2))/delt;
                    %     gD(:,i) =  gD(:,i) + obj.terms{k}.propensity.ODEstoichVector*wp(k);
                    % end
                % end
                A = [gA,gB;gC,gD];
                ANaN = isnan(A);
                if sum(ANaN(:))>0
                    warning('NaNs detected and set to zero. Errors may be present.')
                    A(ANaN) = 0;
                end
            end

        end
        function A = hybridMatrix(obj,t,parameters,v2)
            % Form the infinitesimal generator matrix A
            
            if obj.terms{1}.isFactorizable
                wt = obj.terms{1}.propensity.hybridFactorVector(t,parameters,v2');
                % A = obj.terms{1}.propensity.hybridFactor(t,v2)*obj.terms{1}.matrix;
                A = wt(1)*obj.terms{1}.matrix;
            else
                A = ssit.FspMatrixTerm.generateHybridgMatrixTerm(t, obj.terms{1}.propensity, obj.terms{1}.matrix, parameters, obj.terms{1}.numConstraints, v2);
            end

            for i = 2:length(obj.terms)
                if obj.terms{i}.isFactorizable
                    % A = A + obj.terms{i}.propensity.hybridFactor(t,v2)*obj.terms{i}.matrix;
                    A = A + wt(i)*obj.terms{i}.matrix;
                else
                    A = A + ssit.FspMatrixTerm.generateHybridgMatrixTerm(t, obj.terms{i}.propensity, obj.terms{i}.matrix, parameters, obj.terms{i}.numConstraints, v2);
                end
            end

        end
    end
end
