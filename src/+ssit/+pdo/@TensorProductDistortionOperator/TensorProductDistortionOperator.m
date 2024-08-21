classdef TensorProductDistortionOperator < ssit.pdo.AbstractDistortionOperator
    % (Truncated) Distortion Operator in which the observation channels are
    % are mutually independent discrete random variables and each channel
    % depends on only a single species.
    %
    % Mathematically, this operator can be factorized using the tensor
    % product as
    %   C = C₁ ⊗ C₂ ⊗ ... ⊗ Cₙ
    % where n = number of species = number of observation channels.
    %
    % Parameters
    % ----------
    %
    % conditionalPmfs: cell array of function handles
    %   Conditional probability mass functions. Each of
    %   these functions must take two arguments (x,y) where x is
    %   the true count and y the observed count, and return the
    %   probability of observing y given x. The length of
    %   `conditionalPmfs` must equal the number of species in the
    %   CME model. conditionalPmfs{i} must correspond to the conditional
    %   distribution of observation channel `i` given species `i`.
    %
    % observationDomains: cell array of arrays
    %   Sets of values for observation channels by which to evaluate
    %   probabilities at.
    properties
        conditionalPmfs
        dCdLam
%         observationDomains
    end

    methods
        function obj = TensorProductDistortionOperator(conditionalPmfs,dCdLam)
            %   `obj = TensorProductDistortionOperator(conditionalPmfs)`
            %   constructs a tensor-product distortion operator from a cell
            %   array of conditional probability mass functions `conditionalPmfs` and
            %   single-observation domains `observationDomains`. Each of
            %   these functions must take two arguments `(x,y)` where x is
            %   the true count and y the observed count, and return the
            %   probability of observing y given x. The length of
            %   `conditionalPmfs` must equal the number of species in the
            %   CME model. `observationDomains` must be a cell array where
            %   `observationDomains{i}` is the array of values for the `i`-th
            %   observation channel by which to evaluate probabilities at.
            obj.conditionalPmfs = conditionalPmfs;
            obj.dCdLam = dCdLam;
%             obj.observationDomains = observationDomains;
        end

        function py = computeObservationDist(obj, px, indsIgnore)
            % Compute the probability distribution of
            %distorted single-cell observations using the FSP-approximated
            %probability distribution of true single-cell molecular counts.
            %
            % Parameters
            % ----------
            % px: :mat:class:`~+ssit.@FspVector.FspVector`
            %   Probability distribution of true CME states.
            %
            % Returns
            % -------
            % py: :mat:class:`~+ssit.@FspVector.FspVector`
            % probability distribution of distorted measurements.
            arguments
                obj
                px
                indsIgnore=[] % indices to ignore due to missing observations
            end

            speciesBounds = size(px.data);
            speciesCount = length(speciesBounds);
            pdoFactors = cell(speciesCount, 1);
            kSpecies = 0;

            allPDOsProvided = speciesCount == length(obj.conditionalPmfs);
            % Check to see if all PDOs are provided, or just an orderred
            % subset.

            for iSpecies = 1:speciesCount
                if min(abs(indsIgnore-iSpecies))==0
                    pdoFactors{iSpecies} = ones(1,speciesBounds(iSpecies));
                else
                    if allPDOsProvided
                        kSpecies = iSpecies; % All PDOS are provided.
                    else
                        kSpecies = kSpecies+1; % Use th next in the provided list.
                    end

                    % speciesBounds = size(obj.conditionalPmfs{iSpecies},2)
                    pdoFactors{iSpecies} = obj.conditionalPmfs{kSpecies}(:,1:speciesBounds(iSpecies));
                    nonZeroRows = find(sum(pdoFactors{iSpecies},2)~=0,1,'last');
                    pdoFactors{iSpecies} = obj.conditionalPmfs{kSpecies}(1:nonZeroRows,1:speciesBounds(iSpecies));            
                end
            end
            if numel(px.data)>1
                py = ssit.FspVector(ttm(px.data, pdoFactors));
            else
                py = ssit.FspVector(ttm(full(px.data), pdoFactors));
            end
            py.data = sptensor(py.data);
        end

        function sy = computeObservationDistDiff(obj, px, sx, parameter_idx)
            % Compute the partial derivative of the observation probability distribution
            % with respect to the parameter with index `prameter_idx`, given FSP-approximated probability
            % distribution `px` of true single-cell states and its partial derivative `sx` computed by Forward Sensitivity
            % FSP.
            sy = computeObservationDist(obj, sx);
        end

        function dCdLtimesPx = computeDiffPdoPx(obj, px, ~, parameter_idx)
            % Compute the partial derivative PDO times px.
            speciesBounds = size(px.data);
            speciesCount = length(speciesBounds);
            
            pdoFactors = obj.dCdLam(:,parameter_idx);

            for iSpecies = 1:speciesCount
                if ~isempty(pdoFactors{iSpecies})
                    pdoFactors{iSpecies} = pdoFactors{iSpecies}(:,1:speciesBounds(iSpecies));
                end
            end
            
            pdoParCount = size(pdoFactors,1);
            for iSpecies = 1:pdoParCount
                sz(iSpecies) = ~isempty(pdoFactors{iSpecies});
            end
            idx = find(sz);
            if isempty(idx)
                dCdLtimesPx = [];
            else
                if numel(px.data)>1
                    dCdLtimesPx = ssit.FspVector(ttm(px.data, pdoFactors, idx));
                else
                    dCdLtimesPx = ssit.FspVector(ttm(full(px.data), pdoFactors, idx));
                end
                dCdLtimesPx.data = sptensor(dCdLtimesPx.data);
            end
        end
    end
end


