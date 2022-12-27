function [solutions, constraintBoundsFinal,stateSpace] = adaptiveFspSolve(...
    outputTimes, ...
    initStates,...
    initProbs,...
    stoichMatrix, propensities, ...
    fspTol,...
    constraintFunctions, initialConstraintBounds,...
    verbose,...
    relTol,absTol,odeSolver,...
    stateSpace,usePiecewiseFSP,initApproxSS)
% Approximate the transient solution of the chemical master equation using
% an adaptively expanding finite state projection (FSP).
%
% Parameters
% ----------
%
%   outputTimes: 1-D array
%        time points at which to output the solutions.
%
%   initStates: 2-D array
%       initial states, with shape (number_of_species)x(number_of_states)
%
%   initProbs: 1-D column vector
%       initial probabilities, ordered according to the order of the states in initStates.
%
%   stoichMatrix: 2-D array
%       stoichiometry matrix, of shape (number of species)x(number of reactions)
%
%   propensities: cell of :mat:class:`~+ssit.@Propensity.Propensity` objects
%       the j-th cell entry is the propensity function corresponding to the j-th reaction.
%
%   fspTol: double
%       FSP tolerance.
%
%   constraintFunctions: function handle
%       function to compute the FSP shape constraints. It must be callable with the
%       syntax ``Y = fConstraints(X)`` where X is the array of states arranged
%       column-wise, and each row of the output array Y is a component of the
%       constraint. For example, ``fConstraints = @(X) [X(1,:); X(2,:)]``
%       is a valid constraint function.
%
%   initialConstraintBounds: 1-d column vector
%        initial constraint bounds.
%
%    verbose: logical
%       whether to output information about the intermediate
%       steps to the terminal.
%
%   relTol: double
%       relative tolerance for the ODE solver.
%
%   absTol: double
%    absolute tolerance for the ODE solver.
%
%   odeSolver: string
%    choice of ODE solver for the truncated problem. Must be one of the three
%    options "sundials", "expokit", "matlab", "auto". Default: "auto".
%
%   stateSpace: ssit.FiniteStateSet  
%    State space for FSP projection 
%    (using class ssit.FiniteStateSet)
%
%   usePiecewiseFSP: logical(false) 
%    Option to use piecewise constant FSP for time
%    varying problems (allows use of expokit).
%
%   initApproxSS: logical(false)
%    option to use reflecting boundary condition to
%    estimate the steady state distribution and use that as the initial
%    condiiton.
%
% Returns
% -------
%
%   solutions: cell
%       the i-th entry is a struct that stores the solution at time
%       ``tSpan(i)``. This struct consists of three properties:
%
%           ``states`` -- array of states, arranged column-wise.
%
%           ``p`` -- probability values, with ``p(j) ==`` probability at
%
%           ``states(:,j)``.
%
%           ``sinks`` -- probability mass absorbed at the sink states.
%
%   bConstraintsFinal: column vector
%       final constraint bounds used for generating the FSP.
%
% See Also
% --------
% :mat:func:`~+ssit.+fsp.marginals`
% :mat:class:`~+ssit.@Propensity.Propensity`
%
arguments
    outputTimes (:,1) double
    initStates (:,:) {mustBeInteger(initStates)}
    initProbs (:,1) double
    stoichMatrix (:,:) {mustBeInteger(stoichMatrix)}
    propensities (:, 1)
    fspTol (1,1) double
    constraintFunctions function_handle
    initialConstraintBounds (:,1) double    
    verbose logical=false    
    relTol double=1.0e-4
    absTol double=1.0e-8
    odeSolver string="auto"
    stateSpace =[];
    usePiecewiseFSP=false;
    initApproxSS=false;
end

maxOutputTime = max(outputTimes);
outputTimeCount = length(outputTimes);
outputTimes = unique(outputTimes);
if (outputTimeCount ~= length(outputTimes))
    disp('Warning: input tspan contains repeated elements. Converted tspan to vector of unique values');
    outputTimeCount = length(outputTimes);
end
constraintCount = size(initialConstraintBounds, 1);
constraintBoundsFinal = initialConstraintBounds;
solutions = cell(outputTimeCount, 1);

% First guess of FSP should include the initial condition.
constraintBoundsFinal = max(constraintBoundsFinal,constraintFunctions(initStates));

% Set up the initial state subset
if isempty(stateSpace)
    stateSpace = ssit.FiniteStateSet(initStates, stoichMatrix);
    stateSpace = stateSpace.expand(constraintFunctions, constraintBoundsFinal);
end

% Generate the FSP matrix
stateCount = stateSpace.getNumStates();
Afsp = ssit.FspMatrix(propensities, stateSpace, constraintCount);

% Check that the model propensities are time-invariant
isTimeInvariant = true;
for i = 1:length(propensities)
    f = propensities{i};
    if f.isTimeDependent
        isTimeInvariant = false;
        break;
    end
end

% Use Approximate steady state as initial distribution if requested.
if initApproxSS
    jac = Afsp.createSingleMatrix(-1);
    jac =jac(1:end-constraintCount,1:end-constraintCount);
    jac =jac+diag(sum(jac));
    try
        [eigVec,~] = eigs(jac,1,'smallestabs');
    catch
        try
            [eigVec,~] = eigs(jac,0);
        catch
            try 
                eigVec = null(full(jac));
            catch
                disp('Could not find null space. Using uniform.')
                eigVec = ones(size(jac,1),1);
            end
        end
    end
    solVec = [eigVec/sum(eigVec);zeros(constraintCount,1)];
else % otherwise use user supplied IC.
    solVec = zeros(stateCount + constraintCount, 1);
    solVec(1:size(initStates,2)) = initProbs;
end

tInit =min(outputTimes);
tNow = tInit;
iout = find(outputTimes == tNow, 1, 'first');
if (~isempty(iout))
    solutions{iout} = struct(time=0, ...
        p=packFspSolution(stateSpace.states, solVec(1:stateCount)), ...
        sinks=solVec(stateCount+1:end));
else
    iout = 0;
end

% Determine ODE solver
if odeSolver == "auto"
    if isTimeInvariant 
        odeSolver = "expokit";
    elseif usePiecewiseFSP
        odeSolver = "expokitPiecewise";
    else
        if (exist('FspCVodeMex')~=3)
            odeSolver = "matlab";
        else
            odeSolver = "sundials";
        end
    end       
end

while (tNow < maxOutputTime)
    fspErrorCondition = struct('tStart', tNow,...
        'tFinal', maxOutputTime,...
        'tInit', tInit,...
        'fspTol', fspTol,...
        'nSinks', constraintCount);

    tOut = outputTimes(outputTimes >= tNow);
    % Set up ODE solver data structure for the truncated problem
    if (odeSolver=="sundials")
        jac = [];
        matvec = @(t, p) Afsp.multiply(t, p);
        solver = ssit.fsp_ode_solvers.MexSundials();        
    else
        if isTimeInvariant==1 % If no parameters are functions of time.
            jac = Afsp.createSingleMatrix(tNow);
            matvec = @(~,p) jac*p;            
        elseif usePiecewiseFSP
            jac = @(t)Afsp.createSingleMatrix(t);
            matvec = [];            
        else
            % If time varrying
            matvec = @(t, p) Afsp.multiply(t, p);
            jac = @(t,~)Afsp.createSingleMatrix(t);            
        end

        if odeSolver == "expokit"
            solver = ssit.fsp_ode_solvers.Expokit();
        elseif odeSolver == "expokitPiecewise"
            solver = ssit.fsp_ode_solvers.ExpokitPiecewise();
        else
            solver = ssit.fsp_ode_solvers.OdeSuite(relTol, absTol);
        end
    end

    [~, solutionsNow, fspStopStatus] = solver.solve(tNow, tOut,...
        solVec, ...
        matvec, ...
        jac,...
        fspErrorCondition);

    j = 0;
    while (j < length(solutionsNow))
        if (isempty(solutionsNow{j+1}))
            break;
        end
        j = j +1;
        iout = iout+1;
        solutions{iout} = struct(time=outputTimes(iout),...
            p=packFspSolution(stateSpace.states, solutionsNow{j}(1:stateCount)),...
            sinks=solutionsNow{j}(stateCount+1:end));
    end

    if (j > 0)
        tNow = outputTimes(iout);
        solVecOld = solutionsNow{j};
    else
        solVecOld = solVec;
    end

    if (fspStopStatus.i_expand)
        if (verbose)
            fprintf('expand at t= %2.2e from fsp with size %d \n',...
                fspStopStatus.t_break, stateSpace.getNumStates);
        end

        constraintsToRelax = find(fspStopStatus.sinks*fspErrorCondition.nSinks >=...
            fspStopStatus.error_bound);
        constraintBoundsFinal(constraintsToRelax) = 1.2*constraintBoundsFinal(constraintsToRelax);

        if min(constraintsToRelax)<=2
            stateSpace = ssit.FiniteStateSet(initStates, stoichMatrix);
            stateSpace = stateSpace.expand(constraintFunctions, constraintBoundsFinal);
            warning('Regenerate State Space')
        else
            try
                stateSpace = stateSpace.expand(constraintFunctions, constraintBoundsFinal);
            catch
                stateSpace = ssit.FiniteStateSet(initStates, stoichMatrix);
                stateSpace = stateSpace.expand(constraintFunctions, constraintBoundsFinal);
            end
        end

        try
            Afsp = Afsp.regenerate(propensities, stateSpace, constraintCount);
        catch
            stateSpace = ssit.FiniteStateSet(initStates, stoichMatrix);
            stateSpace = stateSpace.expand(constraintFunctions, constraintBoundsFinal);
            Afsp = Afsp.regenerate(propensities, stateSpace, constraintCount);
        end

        stateCountOld = stateCount;
        stateCount = stateSpace.getNumStates;
        solVec = zeros(stateCount + constraintCount, 1);
        solVec(1:stateCountOld) = solVecOld(1:stateCountOld);
        solVec(stateCount +1: end) = solVecOld(stateCountOld+1:end);
    else
        solVec = solVecOld;
    end
end
end

function f = packFspSolution(states, p)
f = ssit.FspVector(states, p);
end

