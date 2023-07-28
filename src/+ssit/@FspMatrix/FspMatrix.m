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
        function obj = FspMatrix(propensities, stateSet, numConstraints, varNames, modRedTransformMatrices)
            arguments
                propensities
                stateSet
                numConstraints
                varNames =[]
                modRedTransformMatrices =[];
            end
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
            % Returns
            % -------
            %
            %   obj: an instance of this class.
            %
            obj = obj.regenerate(propensities, stateSet, numConstraints, varNames, modRedTransformMatrices);
        end

        function obj = regenerate(obj, propensities, stateSet, numConstraints, varNames, modRedTransformMatrices)
            arguments
                obj
                propensities
                stateSet
                numConstraints
                varNames = []
                modRedTransformMatrices = []
            end
            n_reactions = length(propensities);
            obj.terms = cell(n_reactions, 1);
            for i = 1:n_reactions
                obj.terms{i} = ssit.FspMatrixTerm(propensities{i}, stateSet, numConstraints, varNames, modRedTransformMatrices);
            end
        end

        function w = multiply(obj,t,v)
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

            w = obj.terms{1}.multiply(t,v);
            for i = 2:length(obj.terms)
                w = w + obj.terms{i}.multiply(t,v);
            end
        end

        function A = createSingleMatrix(obj, t)
            % Generate a single MATLAB sparse matrix from the current time
            % `t`. This matrix's left multiplication on a vector `v` of
            % appropriate size will return the same output as when calling
            % `multiply(obj,t, v)`.
            if obj.terms{1}.isFactorizable
                if (obj.terms{1}.isTimeDependent)
                    A = obj.terms{1}.propensity.timeDependentFactor(t)*obj.terms{1}.matrix;
                else
                    A = obj.terms{1}.matrix;
                end
            else
                A = ssit.FspMatrixTerm.generateTimeVaryingMatrixTerm(t, obj.terms{1}.propensity, obj.terms{1}.matrix, obj.terms{1}.numConstraints);
            end

            for i = 2:length(obj.terms)
                if obj.terms{i}.isFactorizable
                    if obj.terms{i}.isTimeDependent
                        A = A + obj.terms{i}.propensity.timeDependentFactor(t)*obj.terms{i}.matrix;
                    else
                        A = A + obj.terms{i}.matrix;
                    end
                else
                    A = A + ssit.FspMatrixTerm.generateTimeVaryingMatrixTerm(t, obj.terms{i}.propensity, obj.terms{i}.matrix, obj.terms{i}.numConstraints);
                end
            end
            
            ANaN = isnan(A);
            if sum(ANaN(:))>0
                warning('NaNs detected and set to zero. Errors may be present.')
                A(ANaN) = 0;
            end
            
        end
    end
end
