classdef (Abstract) AbstractDistortionOperator    
    % This class defines the interface that all concrete probabilistic
    % distortion operators must implement. 
    methods (Abstract)
        % All concrete subclass of AbstractDistortionOperator must
        % implement at least the following three methods

        py = computeObservationDist(obj,px)
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

        sy = computeObservationDistDiff(obj, px, sx, parameter_idx)
        % Compute the partial derivative of the observation probability distribution 
        % with respect to the parameter with index `prameter_idx`, given FSP-approximated probability
        % distribution of true single-cell states and its partial derivative computed by Forward Sensitivity 
        % FSP.
    end
end

