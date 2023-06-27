function [outputs,constraintBounds] = adaptiveFspFiniteDiff(...
    srnModel,...
    parameters,...
    perturbationFactor,...
    outputTimes,...
    initialStates,...
    initialProbabilities,...
    fspTol,...
    constraintFunctions, constraintBounds,...
    verbose,usePiecewiseFSP,initApproxSS,stateSpace,varNames,useParallel,...
    fspSoln)
arguments
    srnModel
    parameters
    perturbationFactor
    outputTimes
    initialStates
    initialProbabilities
    fspTol
    constraintFunctions
    constraintBounds
    verbose=false;
    usePiecewiseFSP=false;
    initApproxSS=false;
    stateSpace=[];
    varNames=[];
    useParallel = false;
    fspSoln = [];
end
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
if ~isempty(varNames)
    propensities = srnModel.createPropensities(parameters,varNames);
else
    propensities = srnModel.createPropensities(parameters);
end

% Find solution to FSP
if isempty(fspSoln)
    [solutions,constraintBounds,stateSpace] = ssit.fsp.adaptiveFspSolve(outputTimes, initialStates,...
        initialProbabilities, stoichMatrix, propensities, fspTol, constraintFunctions, constraintBounds,...
        verbose,1.0e-4,1.0e-8,'auto',[],usePiecewiseFSP,initApproxSS,varNames);
else
    solutions = fspSoln.fsp;
    stateSpace = fspSoln.stateSpace;
    fspTol = inf;
end


dsolutions = cell(length(outputTimes),size(parameters,1));
if useParallel
    parfor i = 1:size(parameters,1)
        dsolutions(:,i) = computeSensitivityForParameter(parameters,perturbationFactor,...
            srnModel,i,constraintBounds,outputTimes,initialStates,initialProbabilities,...
            stoichMatrix,fspTol,constraintFunctions,verbose,...
            stateSpace,usePiecewiseFSP,initApproxSS,varNames);
    end
else
    for i = 1:size(parameters,1)
        dsolutions(:,i) = computeSensitivityForParameter(parameters,perturbationFactor,...
            srnModel,i,constraintBounds,outputTimes,initialStates,initialProbabilities,...
            stoichMatrix,fspTol,constraintFunctions,verbose,...
            stateSpace,usePiecewiseFSP,initApproxSS,varNames);
    end
end

dsolutionsPack = cell(length(outputTimes), 1);
for i = 1:length(outputTimes)
    dsolutionsPack{i} = ssit.FspVector.empty(0);
end
for i = 1:size(parameters,1)
    for j = 1:length(outputTimes)
        try
            dsolutionsPack{j}(i) = dsolutions{j,i};
        catch
            warning('sparse vectors are not the same size')
        end
    end
end

outputs = cell(length(outputTimes),1);
for i = 1:length(outputTimes)
    outputs{i} = struct(...
        time=solutions{i}.time,...
        p=solutions{i}.p,...
        sinks = solutions{i}.sinks,...
        S = dsolutionsPack(j)...
        );
end
end

function dsolutions = computeSensitivityForParameter(parameters,perturbationFactor,...
    srnModel,i,constraintBounds,outputTimes,initialStates,initialProbabilities,...
    stoichMatrix,fspTol,constraintFunctions,verbose,...
    stateSpace,usePiecewiseFSP,initApproxSS,varNames)

dsolutions = cell(length(outputTimes), 1);
for j = 1:length(outputTimes)
    dsolutions{j} = ssit.FspVector.empty(0);
end
parPlus = parameters;
parMinus = parameters;
% Create upper and lower bound for parameters
parPlus{i, 2} = parPlus{i, 2} + perturbationFactor*parPlus{i, 2};
parMinus{i, 2} = parMinus{i, 2} - perturbationFactor*parMinus{i, 2};

propensitiesPlus = srnModel.createPropensities(parPlus);
propensitiesMinus = srnModel.createPropensities(parMinus);

cstBnd2 = 0*constraintBounds;
while max(abs(cstBnd2-constraintBounds))~=0
    %Find FSP solution at upper bound
    cstBnd2 = constraintBounds;
    [solsPlus,constraintBounds,stateSpace] = ssit.fsp.adaptiveFspSolve(outputTimes, initialStates,...
        initialProbabilities, stoichMatrix,...
        propensitiesPlus, fspTol, ...
        constraintFunctions, constraintBounds,...
        verbose,1.0e-4,1.0e-8,'auto',stateSpace,usePiecewiseFSP,initApproxSS,varNames);

    %Find FSP solution at lower bound
    [solsMinus,constraintBounds,stateSpace] = ssit.fsp.adaptiveFspSolve(outputTimes, initialStates,...
        initialProbabilities, stoichMatrix,...
        propensitiesMinus, fspTol,...
        constraintFunctions, constraintBounds,...
        verbose,1.0e-4,1.0e-8,'auto',stateSpace,usePiecewiseFSP,initApproxSS,varNames);
end

% Find the average of the upper and lower bound vectors
for j = 1:length(outputTimes)
    try
        dsolutions{j} = ssit.FspVector((1.0/abs(2*perturbationFactor*parameters{i, 2}))*(solsPlus{j}.p.data - solsMinus{j}.p.data));
    catch
        warning('sparse vectors are not the same size')
    end
end
end
