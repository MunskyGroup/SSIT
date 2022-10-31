classdef (HandleCompatible) FspVector
    % This class represents the truncated probability distributions arising
    % from the FSP-based numerical solution of the CME. Mathematically,
    % such a solution would be a discrete multivariate probability distribution on the
    % N-dimensional integer lattice with a sparse support. Thus, we can
    % think of it in computational term as a sparse N-dimensional array.
    %
    % Parameters
    % ----------
    %
    % dim: integer
    %   Number of dimensions.
    %
    % data: integer
    %   The sparse tensor data.
    %
    % Notes
    % -----
    %
    % In the current version, this object is
    % a thin wrapper for Sandia Tensor Toolbox [#f1]_ `sptensor` object to support
    % zero-based multi-indices.
    %
    % .. rubric:: Footnotes
    %
    % .. [#f1] https://gitlab.com/tensors/tensor_toolbox    
    properties
        dim
        data
    end

    methods
        function obj = FspVector(varargin)
            % Construct an instance of `FspVector`.
            %
            %   FspVector(DIM) creates an empty object for storing
            %   values of a discrete multivariate probability distribution with DIM variables.
            %
            %   FspVector(T) where T is a `sptensor` object
            %   creates an `FspVector` wrapper for T that supports
            %   zero-based indexing.
            %
            %   FspVector(X, Y) creates an object with variable values X
            %   and function values Y. Y must be a column vector. Each
            %   column of X is a multi-dimensional vector of non-negative
            %   integers and the corresponding entry in Y is the function
            %   value.

            if nargin == 1
                if isa(varargin{1}, 'sptensor') || isa(varargin{1}, 'tensor')
                    obj.data = varargin{1};
                    obj.dim = length(size(obj.data));
                else
                    obj.dim = varargin{1};
                    obj.data = sptensor();
                end
            else
                X = varargin{1};
                Y = varargin{2};
                obj.dim = size(X, 1);
                obj.data = sptensor();
                obj = obj.setValues(X, Y);
            end
        end

        function fout = sumOver(obj, axes)
            %SUMOVER Constructing a new `FspVector` with fewer variables
            % by summing over the axes specified in the input.
            newData = squeeze(sptensor(collapse(obj.data, axes)));
            fout = ssit.FspVector(newData);
        end

        function y = sum(obj)
            %SUM returns the sum of all entries.
            y = squeeze(collapse(obj.data, 1:obj.dim));
        end

        function obj = setValues(obj, Xs, fXs)
            % Set the values of the `FspVector` at the specified locations
            % `Xs` with new values given in `fXs`. `Xs` must be aranged
            % column-wise so that `obj.dim == size(Xs, 1)` and `fXs` must
            % be a column vector.
            if obj.dim ~= size(Xs, 1)
                error("Dimensionality of input data incompatible with object.\n");
            end
            inputIdxs = Xs'+1;
            try
                obj.data(inputIdxs) = fXs;
            catch
                warning('NaN error encountered - setting solution to zero')
                obj.data(inputIdxs) = 0;
            end
        end
    end
end

