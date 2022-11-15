function [outputs,constraintBounds] = adaptiveFspFiniteDiff(...
    srnModel,...
    parameters,...
    perturbationFactor,...
    outputTimes,...
    initialStates,...
    initialProbabilities,...
    fspTol,...
    constraintFunctions, constraintBounds,...
    verbose)
%
% Approximate the partial derivatives of the transient solution of the 
% chemical master equation using finite difference. 
%
% Parameters
% ----------
%
%   srnModel: :mat:class:`~+ssit.@SrnModel.SrnModel`
%       Stochastic reaction network model.
%       
%   parameters: containers.Map
%       Dictionary containing parameter-value pairs.
%
%   outputTimes: vector of scalars
%       Array of output time points.
%
%   initial_states: 2-D array
%       Array of initial states, each state is a column vector.
%
%   initial_prob: column vector of scalars
%       Initial probabilities.
%
%   fspTol: scalar 
%       FSP error tolerance.
%
%   constraintFunctions: function handle 
%       Function to compute the left-hand-side of FSP shape constraints.
%
%   constraintBounds: 1-d column vector 
%       Initial constraint bounds.
%
%   verbose: true/false 
%       Whether to output information about the intermediate steps to the command window.
%
% Returns
% -------
%
%   outputs: cell array
%       Each element `outputs{j}` corresponds to the solution at
%       `outputTimes{j}` and is a struct with fields:
%       - time: time of the solution.
%       - p: probability distribution stored as a :mat:class:`~+ssit.@FspVector.FspVector` object.
%       - sinks: sink state probabilitities.
%       - S: array of :mat:class:`~+ssit.@FspVector.FspVector` objects,
%       each is the partial derivative of `p` with respect to a model
%       parameter.



% Obtain Stoichiometry from model
stoichMatrix = srnModel.stoichiometry;

% Obtain CME solution at paramter_values
propensities = srnModel.createPropensities(parameters);

% Find solution to FSP
[solutions,constraintBounds] = ssit.fsp.adaptiveFspSolve(outputTimes, initialStates,...
    initialProbabilities, stoichMatrix, propensities, fspTol, constraintFunctions, constraintBounds,...
    verbose);

dsolutions = cell(length(outputTimes), 1);
for i = 1:length(outputTimes)
    dsolutions{i} = ssit.FspVector.empty(0);
end
for i = 1:size(parameters,1)
    parPlus = parameters;
    parMinus = parameters;
    % Create upper and lower bound for parameters
    parPlus{i, 2} = parPlus{i, 2} + perturbationFactor*parPlus{i, 2};
    parMinus{i, 2} = parMinus{i, 2} - perturbationFactor*parMinus{i, 2};

    propensitiesPlus = srnModel.createPropensities(parPlus);
    propensitiesMinus = srnModel.createPropensities(parMinus);

    %Find FSP solution at upper bound
    [solsPlus,constraintBounds] = ssit.fsp.adaptiveFspSolve(outputTimes, initialStates,...
        initialProbabilities, stoichMatrix,...
        propensitiesPlus, fspTol, ...
        constraintFunctions, constraintBounds,...
        verbose);

    %Find FSP solution at lower bound
    [solsMinus,constraintBounds] = ssit.fsp.adaptiveFspSolve(outputTimes, initialStates,...
        initialProbabilities, stoichMatrix,...
        propensitiesMinus, fspTol,...
        constraintFunctions, constraintBounds,...
        verbose);

    % Find the average of the upper and lower bound vectors
    for j = 1:length(outputTimes)
        try
            dsolutions{j}(i) = ssit.FspVector((1.0/abs(2*perturbationFactor*parameters{i, 2}))*(solsPlus{j}.p.data - solsMinus{j}.p.data));
        catch
            1+1
        end
    end
end
outputs = cell(length(outputTimes),1);
for i = 1:length(outputTimes)
    outputs{i} = struct(...
        time=solutions{i}.time,...
        p=solutions{i}.p,...
        sinks = solutions{i}.sinks,...
        S = dsolutions{i}...
        );
end
end

