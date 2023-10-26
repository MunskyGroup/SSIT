function [outputs, constraintBounds, stateSpace] = adaptiveFspForwardSens(outputTimes,...
    initialStates,...
    initialProbabilities, initialSensitivities,...
    stoichMatrix, ...
    propensities, parameters, propensityDerivatives, computableSensitivities,...
    fspTol,...
    varNames, modRedTransformMatrices, ...
    constraintFunctions, initialConstraintBounds,...
    verbose, useMex,...
    relTol,...
    absTol, ...
    stateSpace)
arguments
    outputTimes
    initialStates
    initialProbabilities
    initialSensitivities
    stoichMatrix
    propensities
    parameters
    propensityDerivatives
    computableSensitivities
    fspTol
    varNames
    modRedTransformMatrices
    constraintFunctions
    initialConstraintBounds
    verbose
    useMex
    relTol
    absTol
    stateSpace =[];
end
% Compute and outputs the solution and sensitivitiy vectors of the CME at the user-input timepoints.
%
% Parameters
% ----------
%
%   outputTimes: array of scalars
%       output time points.
%
%   initialStates: 2-D array
%       array of initial states, each state is a column vector.
%
%   initialProbabilities: vector of scalars
%       column vector of initial probabilities.
%
%   initialSensitivities:
%       initial sensitivities of the initial states, arranged as
%       a long column vector with length n_init * n_par where n_init is the
%       length of initialProbabilities and n_par the number of parameters.
%
%   stoichMatrix: matrix
%       Net stoichiometry matrix.
%
%   propensities: 1-D cell array
%       cell of instances of :mat:class:`~+ssit.@Propensity.Propensity`.
%
%   propensityDerivatives: 2-D cell array
%       cell of propensity objects, but for the partial
%       derivatives of the propensities with respect to parameters.
%       propensityDerivatives{i, j} stores the derivative of the i-th propensity wrt to
%       the j-th parameter.
%
%   fspTol: FSP error tolerance.
%
%   constraintFunctions: function handle
%       function to compute the FSP shape constraints. It must be callable with the
%       syntax ``Y = fConstraints(X)`` where X is the array of states arranged
%       column-wise, and each row of the output array Y is a component of the
%       constraint. For example, ``fConstraints = @(X) [X(1,:); X(2,:)]``
%       is a valid constraint function.
%
%   constraintBounds: column vector of scalars
%        initial constraint bounds.
%
%   verbose: true/false
%       whether to output intermediate status.
%
%   useMex: (true/false)
%       whether to use the MEX Forward Sensitivity solver.

if (exist('FspSensCVodeMex') ~=3)
    useMex = false;
end

% Check that the model propensities are time-invariant
isTimeInvariant = true;
for i = 1:length(propensities)
    f = propensities{i};
    if f.isTimeDependent
        isTimeInvariant = false;
        break;
    end
end

% check input sizes
parameterCount = sum(computableSensitivities);

% if (reactionCount ~= size(propensityDerivatives,1))
%     error('Number of rows of propensityDerivatives cell must be the same as number of propensities');
% end
if (size(initialProbabilities, 1)*parameterCount ~= size(initialSensitivities,1))
    error('Input probabilities and sensitivities must have compatible size.');
end

tFinal = max(outputTimes);
tOutputCount = length(outputTimes);
outputTimes = unique(outputTimes);
if (tOutputCount ~= length(outputTimes))
    disp('Warning: input tspan contains repeated elements. Converted tspan to vector of unique values');
    tOutputCount = length(outputTimes);
end
constraintCount = size(initialConstraintBounds, 1);
constraintBounds = initialConstraintBounds;
outputs = cell(tOutputCount, 1);

% Set up the initial state subset
if isempty(stateSpace)
    stateSpace = ssit.FiniteStateSet(initialStates, stoichMatrix);
    stateSpace = stateSpace.expand(constraintFunctions, constraintBounds);
end
stateCount = stateSpace.getNumStates();

% Generate the time-varying FSP operator
fspMatrix = ssit.FspMatrix(propensities, [parameters{:,2}]', stateSpace, constraintCount, varNames, modRedTransformMatrices, true);

probabilityVec = zeros(stateCount + constraintCount, 1);
probabilityVec(1:size(initialStates,2)) = initialProbabilities;
sensitivityVecs = zeros((stateCount + constraintCount)*parameterCount, 1);
for i = 1:parameterCount
    sensitivityVecs(1 + (i-1)*length(initialProbabilities):i*length(initialProbabilities)) ...
        = initialSensitivities(1 + (i-1)*length(initialProbabilities):i*length(initialProbabilities));
end
sensitivityVecs = reshape(sensitivityVecs, stateCount + constraintCount, parameterCount);

tInit =min(outputTimes);
tNow = tInit;
iout = find(outputTimes == tNow, 1, 'first');
if (~isempty(iout))
    outputs{iout} = packForwardSensFspSolution(0.0,...
        stateSpace.states,...
        probabilityVec, sensitivityVecs);
else
    iout = 0;
end

while (tNow < tFinal)
    stop_cond = struct('t_start', tNow, 't_final', tFinal,...
        'fspTol', fspTol, 'n_sinks', constraintCount, 'n_states', stateCount, 'tInit', tInit);
    % Solve the forward sensitivity FSP ODE system
    if (useMex)
        matvec = @(t, p) fspMatrix.multiply(t, p);
        dmatvec = cell(parameterCount, 1);
        for i = 1:parameterCount
            dmatvec{i} = @(t, v) fspMatrixDiff{i}.multiply(t, v);
        end

        f_out = @(t, p, dp) struct('p', p, 'dp', {dp});

        [outputs_current, stopStatus] = FspSensCVodeMex(tNow, outputTimes(outputTimes>tNow), matvec, dmatvec, ...
            probabilityVec, sensitivityVecs, f_out, @fspErrorSundialsEvent, stop_cond);
        j = 0;
        while (j < length(outputs_current))
            if (isempty(outputs_current{j+1}))
                break;
            end
            j = j +1;
            iout = iout+1;
            DP = cell2mat(outputs_current{j}.dp');
            outputs{iout} = packForwardSensFspSolution(outputTimes(iout),...
                stateSpace.states,...
                outputs_current{j}.p, DP);
        end
        
    else
        y0 = zeros(length(probabilityVec)*(parameterCount+1), 1);
        y0(1:stateCount+constraintCount) = probabilityVec;
        for j = 1:parameterCount
            y0(j*(stateCount+constraintCount)+1:(j+1)*(stateCount+constraintCount)) = ...
                sensitivityVecs(:,j);
        end
        
        if isTimeInvariant
            % Expokit Solution
            jac = forwardSensJac(tNow, y0, fspMatrix, [parameters{:,2}]');
            tryAgain=1;
            m=30;
            indSinks =[stop_cond.n_states+1:stop_cond.n_states+stop_cond.n_sinks];
            while tryAgain==1
                [~, ~, ~, tout, outputs_current, ~, tryAgain, te, ye] =...
                    ssit.fsp_ode_solvers.expv_modified(outputTimes(end), jac,...
                    y0, 1e-8, m, [], outputTimes(outputTimes>=tNow), fspTol, indSinks, tNow, stop_cond);
                if tryAgain==0;break;end
                if m>300
                    warning('Expokit expansion truncated at 300');
                    [~, ~, ~, tout, outputs_current, ~, ~, te, ye] =...
                        ssit.fsp_ode_solvers.expv_modified(outputTimes(end), jac, y0, fspTol/1e5, m,...
                        [], outputTimes(outputTimes>=tNow), fspTol, indSinks, tNow);
                    tryAgain=0;
                    break;
                end
                m=m+5;
            end
            %%
        else
            % ODE Solver
            ode_event = @(t,ps) fspErrorEvent(t, ps, stop_cond);
            jac = @(t,ps,dpsdt) forwardSensJac(t, ps, fspMatrix, fspMatrixDiff(indsCompSens));
            ode_rhs = @(t, ps) forwardSensRHS(t, ps, fspMatrix, fspMatrixDiff(indsCompSens), stateCount, constraintCount, parameterCount);
            ode_opts = odeset(Events=ode_event, Jacobian=jac, RelTol=relTol, AbsTol=absTol);
            [tout, outputs_current, te, ye, ~] = ...
                ode23s(ode_rhs, outputTimes(outputTimes>=tNow), y0, ode_opts);
        end

        if length(tout)<2||(~isempty(te)&&te<tout(2))
            outputs_current = outputs_current(1);
            tout = tout(1);
        end

        ikeep = ismember(tout, outputTimes) & (tout > tNow);
        outputs_current = outputs_current(ikeep, :)';

        if (~isempty(te))
            sinks = ye(stateCount+1:stateCount+constraintCount);
            errorBound = stop_cond.fspTol*(te-stop_cond.t_start)/(stop_cond.t_final-stop_cond.t_start);
            stopStatus = struct(iExpand=true, tBreak=te, sinks=sinks, errorBound=errorBound);
        else
            stopStatus = struct(iExpand=false, tBreak=[], sinks=[], errorBound=[]);
        end
        % Enter intermediate solutions to the Outputs container
        j = 0;
        while (j < size(outputs_current, 2))
            j = j+1;
            iout = iout+1;
            DP = reshape(outputs_current(stateCount+constraintCount+1:end,j), stateCount+constraintCount, parameterCount);
            outputs{iout} = packForwardSensFspSolution(outputTimes(iout),...
                stateSpace.states,...
                outputs_current(1:stateCount+constraintCount,j), DP);
        end
    end

    if (j > 0)
        tNow = outputTimes(iout);
        solVecOld = outputs_current(1:stateCount+constraintCount,j);
        sensitivityVecsOld = reshape(outputs_current(stateCount+constraintCount+1:end,j), stateCount+constraintCount, parameterCount);
    else
        solVecOld = probabilityVec;
        sensitivityVecsOld = sensitivityVecs;
    end

    if (stopStatus.iExpand)
        if (verbose)
            fprintf('expand at t= %2.2e from fsp with size %d \n', stopStatus.tBreak, stateSpace.getNumStates);
        end

        i_relax = find(stopStatus.sinks*constraintCount >= 0.8*stopStatus.errorBound);
        constraintBounds(i_relax) = 1.2*constraintBounds(i_relax);

        stateSpace = stateSpace.expand(constraintFunctions, constraintBounds);

        fspMatrix = ssit.FspMatrix(propensities, [parameters{:,2}]', stateSpace, constraintCount, varNames, modRedTransformMatrices, true);

        stateCountOld = stateCount;
        stateCount = stateSpace.getNumStates;

        probabilityVec = zeros(stateCount + constraintCount, 1);
        sensitivityVecs = zeros((stateCount + constraintCount), parameterCount);

        probabilityVec(1:stateCountOld) = solVecOld(1:stateCountOld);
        probabilityVec(stateCount +1: end) = solVecOld(stateCountOld+1:end);
        for i = 1:parameterCount
            sensitivityVecs(1:stateCountOld, i) ...
                = sensitivityVecsOld(1:stateCountOld, i);
        end
    end
end
end


function [stopHere, stopStatus] = fspErrorSundialsEvent(t, p, stop_cond)
% This function allows the Sundials MEX solver to detect when the FSP error exceeds tolerance.
%
% Parameters
% ----------
%
%   t: scalar
%     current time.
%
%   p: vector of real scalars
%     current solution vector (including sink state probabilities at the end).
%
%   stop_cond: struct
%      data structure for FSP error tolerance. The fields are `fspTol`,
%      `t_final`, `n_sinks` for the FSP tolerance, final time, number of
%      sinks.
%
% Returns
% -------
%
%   stopHere: scalar
%       is 0 if FSP error is still acceptable, otherwise 1.
%
%   stopStatus: struct
%       FSP error status at time `t`. Fields are
%           `iExpand`: true if the FSP needs expansion, false otherwise.
%           `tBreak`: time the FSP needs to break.
%           `sinks`: probabilities of sink states.
%           `errorBound`: The error bound condition that the solution needs
%           to satisfy.
%
%
sinks = p(end+1-stop_cond.n_sinks:end);
errorBound = stop_cond.fspTol*(t-stop_cond.t_start)/(stop_cond.t_final-stop_cond.t_start);
if (max(sinks)*stop_cond.n_sinks > errorBound)
    stopHere = double(1);
    stopStatus = struct(iExpand=true, tBreak=t, sinks=sinks, errorBound=errorBound);
else
    stopHere = double(0);
    stopStatus = struct(iExpand=false, tBreak=t, sinks=sinks, errorBound=errorBound);
end
end

function [val, terminal, direction] = fspErrorEvent(t, ps, stop_cond)
% This function allows MATLAB's ODE Suite solvers to detect the event when
% the FSP error exceeds tolerance.
%
% Parameters
% ----------
%
%   t: scalar
%     current time.
%
%   ps: vector of real scalars
%     current solution vector (including sink state probabilities at the end).
%
%   stop_cond: struct
%      data structure for FSP error tolerance. The fields are `fspTol`,
%      `t_final`, `n_sinks`, `n_states` for the FSP tolerance, final time, number of
%      sinks, number of FSP states.
%
% Returns
% -------
%
%   val: scalar
%       The value for FSP error subtracted by the required error tolerance.
%       If `val` is negative, the ODE solver continues. Otherwise, it will reduce
%       stepsize until `val` equals zero and returns.
%
%   terminal: scalar
%       Always returns 1. This value tells MATLAB's solver that it must
%       terminate when FSP error reaches tolerance.
%
%   direction: scalar
%       Always returns 1. This tells MATLAB that the FSP error will reach
%       the tolerance from below.
%
% See also
% ---------
% Find more details in the official documentation of MATLAB's ODE solvers.
% Especially, look up for "event detection".

nstates = stop_cond.n_states;
sinks = ps(nstates+1:nstates+stop_cond.n_sinks);
error_bound = stop_cond.fspTol*(t-stop_cond.t_start)/(stop_cond.t_final-stop_cond.t_start);

val = max(sinks)*stop_cond.n_sinks - error_bound;
terminal = 1;
direction = 1;
end

function y = forwardSensRHS(t, ps, A, dA, n_states, n_sinks, n_pars)
y = zeros((n_states+n_sinks)*(n_pars+1),1);
y(1:n_states+n_sinks) = A.multiply(t, ps(1:n_states+n_sinks));
for j = 1:n_pars
    y(j*(n_states+n_sinks)+1:(j+1)*(n_states+n_sinks)) = ...
        A.multiply(t, ps(j*(n_states+n_sinks)+1:(j+1)*(n_states+n_sinks))) + ...
        dA{j}.multiply(t, ps(1:n_states+n_sinks));
end
end

function J = forwardSensJac(t, x, A, parameters, modRedTransformMatrices)
arguments
    t
    x
    A
    parameters
    modRedTransformMatrices = []
end
[Amerged,AmergedSens] = A.createSingleMatrix(t,parameters,modRedTransformMatrices,true);

n = size(Amerged,1);
npars = size(parameters,1);
nzbound = nnz(Amerged)*(npars+1);
J = spalloc(n, n, nzbound);
J(1:n, 1:n) = Amerged;
for j=1:npars
    J(j*n+1:(j+1)*n, j*n+1:(j+1)*n) = Amerged;
    %BRIAN     J(j*n+1:(j+1)*n, (j-1)*n+1:j*n) = dA{j}.createSingleMatrix(t);
    % J(j*n+1:(j+1)*n, 1:n) = dA{j}.createSingleMatrix(t);
    J(j*n+1:(j+1)*n, 1:n) = AmergedSens{j};
end
end

function f = packForwardSensFspSolution(time, states, pvec, svecs)
stateCount = size(states,2);
sensIndices = ssit.FspVector.empty(0);
for ip = 1:size(svecs, 2)
    sensIndices(ip) = ssit.FspVector(states, svecs(1:stateCount,ip));
end
f = struct( time=time,...
    p=ssit.FspVector(states, pvec(1:stateCount)),...
    sinks=pvec(stateCount+1:end),...
    S=sensIndices,...
    dsinks=svecs(stateCount+1:end,:));
end
