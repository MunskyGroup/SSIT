function [solutions, constraintBoundsFinal, stateSpace] = adaptiveFspSolve(...
    outputTimes, ...
    initStates,...
    initProbs,...
    stoichMatrix, ...
    propensities, ...
    parameters, ...
    fspTol,...
    constraintFunctions, initialConstraintBounds,...
    verbose,...
    relTol,absTol,odeSolver,...
    stateSpace,usePiecewiseFSP,initApproxSS,speciesNames,...
    useReducedModel,modRedTransformMatrices,...
    useHybrid,hybridOptions,...
    fEscape,bEscape,...
    constantJacobian,constantJacobianTime)
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
%   speciesNames: cell vector
%    list of names of species
%
%   useReducedModel: logical(false)
%    option to use a projection-based model reduction scheme
%
%   modRedTransformMatrices: structure(optional)
%    model reduction transformation matrices.
%
%   useHybrid: logical(false)
%    option to use a hybrid model where some upstrame species are treated
%    as ODEs
%
%   hybridOptions: structure
%    structure containing the field 'upstreamODEs' which specifies the
%    names of all species that are treated as upstream ODEs.
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
    initStates (:,:) %{mustBeInteger(initStates)}
    initProbs (:,1) double
    stoichMatrix (:,:) {mustBeInteger(stoichMatrix)}
    propensities (:, 1)
    parameters 
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
    speciesNames=[];
    useReducedModel=false;
    modRedTransformMatrices=[];
    useHybrid = false;
    hybridOptions = [];
    fEscape = []
    bEscape = []
    constantJacobian = false;
    constantJacobianTime = NaN;
end

maxOutputTime = max(outputTimes);
outputTimeCount = length(outputTimes);
outputTimes = unique(outputTimes);
if (outputTimeCount ~= length(outputTimes))
    disp('Warning: input tspan contains repeated elements. Converted tspan to vector of unique values');
    outputTimeCount = length(outputTimes);
end

% If it is a hybrid model, construct the deterministic and stochastic
% stoichiometry vectors.
if useHybrid
    jStochastic = find(~contains(speciesNames,hybridOptions.upstreamODEs));
    jUpstreamODE = find(contains(speciesNames,hybridOptions.upstreamODEs));
    odeStoichs = stoichMatrix(jUpstreamODE,:);
    stoichMatrix = stoichMatrix(jStochastic,:);
    initODEs = initStates(jUpstreamODE,1);
    initStates = initStates(jStochastic,:);
    speciesNames = speciesNames(jStochastic);
    numODEs = length(jUpstreamODE);
else
    numODEs = 0;
end

if ~isempty(fEscape)
    initialConstraintBounds = [initialConstraintBounds;bEscape];
    constraintFunctions = @(x)[constraintFunctions(x);fEscape(x)];
    nEscapeSinks = length(bEscape);
else
    nEscapeSinks = 0;
end
constraintCount = size(initialConstraintBounds, 1);
constraintBoundsFinal = initialConstraintBounds;
    
solutions = cell(outputTimeCount, 1);

% First guess of FSP should include the initial condition.
constraintBoundsFinal = max([constraintBoundsFinal,constraintFunctions(initStates)],[],2);

% Set up the initial state subset
if isempty(stateSpace)
    stateSpace = ssit.FiniteStateSet(initStates, stoichMatrix);
    stateSpace = stateSpace.expand(constraintFunctions, constraintBoundsFinal);
end

% Generate the FSP matrix
stateCount = stateSpace.getNumStates();
if useHybrid
    AfspFull = ssit.FspMatrix(propensities, parameters, stateSpace, constraintCount, speciesNames, modRedTransformMatrices);    
elseif useReducedModel
    AfspRed = ssit.FspMatrix(propensities, parameters, stateSpace, constraintCount, speciesNames, modRedTransformMatrices);
else
    AfspFull = ssit.FspMatrix(propensities, parameters, stateSpace, constraintCount, speciesNames);
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

% Use Approximate steady state as initial distribution if requested.
if initApproxSS
    if useHybrid
        FUN = @(v)odeStoichs*generate_propensity_vector(0, v, zeros(length(jStochastic),1), propensities, parameters);
        OPTIONS = optimoptions('fsolve','display','none',...
            'OptimalityTolerance',1e-8,'MaxIterations',2000);
        x0b = fsolve(FUN,initODEs,OPTIONS);
        FUN = @(t,v)odeStoichs*generate_propensity_vector(0, v, zeros(length(jStochastic),1), propensities, parameters);
        [~,ode_solutions] = ode45(FUN,max(outputTimes)*[0,500,1000],x0b);
        initODEs = ode_solutions(end,:)';

        jac = AfspFull.createJacHybridMatrix(0, initODEs, parameters, length(hybridOptions.upstreamODEs), true);

    else
        if useReducedModel
            AfspFull = ssit.FspMatrix(propensities, parameters, stateSpace, constraintCount, speciesNames);
        end

        jac = AfspFull.createSingleMatrix(outputTimes(1)-1e-6,parameters);
    end
    jac = jac(1:end-constraintCount,1:end-constraintCount);

    try
        warning('off')
        [eigVec,~] = eigs(jac,1,'smallestabs');
    catch
        try
            [eigVec,~] = eigs(jac,1);
        catch
            try
                eigVec = null(full(jac));
            catch
                disp('Could not find null space. Using uniform.')
                eigVec = ones(size(jac,1),1);
            end
        end
    end
    if useReducedModel
        solVec = eigVec/sum(eigVec);
        solVec = modRedTransformMatrices.phi_inv*solVec;
    else
        solVec = [eigVec/sum(eigVec);zeros(constraintCount,1)];
    end
    
else % otherwise use user supplied IC.
    if useReducedModel
        solVec = zeros(stateCount, 1);
        solVec(1:size(initStates,2)) = initProbs;
        solVec = modRedTransformMatrices.phi_inv*solVec;
    else
        solVec = zeros(stateCount + constraintCount, 1);
        solVec(1:size(initStates,2)) = initProbs;
    end
end
if useHybrid
    solVec = [solVec;initODEs];
end

% Write initial condition into results structure
tInit =min(outputTimes);
tNow = tInit;
iout = find(outputTimes == tNow, 1, 'first');
if (~isempty(iout))
    if useReducedModel
        if ~isempty(modRedTransformMatrices.phiPlot)
            p = real(modRedTransformMatrices.phiPlot*solVec);
        else
            p = real(modRedTransformMatrices.phi*solVec);
        end
        p = p/sum(p);
        p = max(p,0);
        solutions{iout} = struct(time=0, ...
            p=packFspSolution(stateSpace.states, p(1:stateCount)), ...
            escapeProbs=[], ...
            sinks=[]);
    elseif useHybrid
        solutions{iout} = struct(time=0, ...
            p=packFspSolution(stateSpace.states, solVec(1:stateCount)), ...
            sinks=solVec(stateCount+1:end-nEscapeSinks-numODEs),...
            escapeProbs=solVec(end-nEscapeSinks-numODEs+1:end-numODEs),...
            upstreamODEs=solVec(end-numODEs+1:end));
    else
        solutions{iout} = struct(time=0, ...
            p=packFspSolution(stateSpace.states, solVec(1:stateCount)), ...
            sinks=solVec(stateCount+1:end-nEscapeSinks), ...
            escapeProbs=solVec(end-nEscapeSinks+1:end));
    end
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

% While loop for solving the FSP over each time step.
while (tNow < maxOutputTime)
    fspErrorCondition = struct('tStart', tNow,...
        'tFinal', maxOutputTime,...
        'tInit', tInit,...
        'fspTol', fspTol,...
        'nSinks', constraintCount, ...
        'nEscapeSinks', nEscapeSinks,...
        'numODEs', numODEs);

    tOut = outputTimes(outputTimes >= tNow);
    % Set up ODE solver data structure for the truncated problem
    if (odeSolver=="sundials")
        error('Sundials not implemented.')
%         jac = [];
%         matvec = @(t, p) Afsp.multiply(t, p);
%         solver = ssit.fsp_ode_solvers.MexSundials();        
    else
        if useReducedModel
            fspErrorCondition.fspTol = inf;
            if isTimeInvariant==1 % If no parameters are functions of time.
                jac = sparse(AfspRed.createSingleMatrix(tNow, parameters,modRedTransformMatrices));
                matvec = @(~,p) jac*p;
            elseif usePiecewiseFSP
                if constantJacobian
                    jac = sparse(AfspRed.createSingleMatrix(constantJacobianTime, parameters,modRedTransformMatrices));
                else
                    jac = @(t)AfspRed.createSingleMatrix(t, parameters,modRedTransformMatrices);
                end
                matvec = [];
            else
                % If time varying
                if constantJacobian
                    jac = sparse(AfspRed.createSingleMatrix(constantJacobianTime, parameters,modRedTransformMatrices));
                    matvec = @(~, p)jac*p;
                else
                    matvec = @(t, p) AfspRed.multiply(t, p, parameters);
                    jac = @(t,~)AfspRed.createSingleMatrix(t, parameters,modRedTransformMatrices);
                end
            end

            if odeSolver == "expokit"
                solver = ssit.fsp_ode_solvers.Expokit(20,1e-8,'expv_modified');
            elseif odeSolver == "expokitPiecewise"
                solver = ssit.fsp_ode_solvers.ExpokitPiecewise();
            else
                solver = ssit.fsp_ode_solvers.OdeSuite(relTol, absTol);
            end
        else
            if isTimeInvariant==1 % If no parameters are functions of time.
                % if ~isempty(parameters)
                    jac = AfspFull.createSingleMatrix(tNow+1e-6, parameters);
                % else
                    % jac = AfspFull.createSingleMatrix(tNow);
                % end
                matvec = @(~,p) jac*p;
            elseif usePiecewiseFSP
                if constantJacobian
                    jac = sparse(AfspFull.createSingleMatrix(constantJacobianTime, parameters));
                    matvec = [];
                else
                    jac = @(t)AfspFull.createSingleMatrix(t, parameters);
                    matvec = [];
                end

            elseif useHybrid
                matvec = @(t, p) full(AfspFull.hybridRHS(t, p, parameters, hybridOptions.upstreamODEs));
                % jac = [];  % Jacobian calculation is not currently available for hybrid models.
                jac = @(t, p)AfspFull.createJacHybridMatrix(t, p, parameters, length(hybridOptions.upstreamODEs));
            else
                if constantJacobian
                    jac = sparse(AfspFull.createSingleMatrix(constantJacobianTime, parameters));
                    matvec = @(~, p)jac*p;
                else
                    matvec = @(t, p) AfspFull.multiply(t, p, parameters);
                    jac = @(t,~)AfspFull.createSingleMatrix(t, parameters);
                end
            end

            if odeSolver == "expokit"
                solver = ssit.fsp_ode_solvers.Expokit(30,absTol);
            elseif odeSolver == "expokitPiecewise"
                if constantJacobian
                    solver = ssit.fsp_ode_solvers.Expokit(30,absTol);
                else
                    solver = ssit.fsp_ode_solvers.ExpokitPiecewise();
                end
            else
%                 if ~isempty(hybridOptions)
%                     solver = ssit.fsp_ode_solvers.OdeSuite(relTol, absTol);
%                 else
                    solver = ssit.fsp_ode_solvers.OdeSuite(relTol, absTol);
%                 end
            end
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
        if useReducedModel
            if ~isempty(modRedTransformMatrices.phiPlot)
                p = real(modRedTransformMatrices.phiPlot*solutionsNow{j});
            else
                p = real(modRedTransformMatrices.phi*solutionsNow{j});
            end
            p = p/sum(p);
            p = max(p,0);
            solutions{iout} = struct(time=outputTimes(iout), ...
                p=packFspSolution(stateSpace.states, p(1:stateCount)), ...
                sinks=[],...
                escapeProbs=[]);
        elseif useHybrid
            solutions{iout} = struct(time=outputTimes(iout),...
                p=packFspSolution(stateSpace.states, solutionsNow{j}(1:stateCount)),...
                sinks=solutionsNow{j}(stateCount+1:end-nEscapeSinks-numODEs), ...
                escapeProbs=solutionsNow{j}(end-nEscapeSinks-numODEs+1:end-numODEs),...
                upstreamODEs=solutionsNow{j}(end-numODEs+1:end));
        else
            solutions{iout} = struct(time=outputTimes(iout),...
                p=packFspSolution(stateSpace.states, solutionsNow{j}(1:stateCount)),...
                escapeProbs=solutionsNow{j}(end-nEscapeSinks+1:end),...
                sinks=solutionsNow{j}(stateCount+1:end-nEscapeSinks));
        end
    
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

        if fspStopStatus.error_bound>0
            constraintsToRelax = find(fspStopStatus.sinks(1:end-fspErrorCondition.nEscapeSinks)*...
                (fspErrorCondition.nSinks-fspErrorCondition.nEscapeSinks) >=...
                fspStopStatus.error_bound(end));
        else
            constraintsToRelax = find(fspStopStatus.sinks(1:end-fspErrorCondition.nEscapeSinks)*...
                (fspErrorCondition.nSinks-fspErrorCondition.nEscapeSinks) > 0);
        end
        constraintBoundsFinal(constraintsToRelax) = 1.2*constraintBoundsFinal(constraintsToRelax);

        if min(constraintsToRelax)<=size(stoichMatrix,1)
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
            AfspFull = AfspFull.regenerate(propensities, parameters, stateSpace, constraintCount,speciesNames);
        catch
            stateSpace = ssit.FiniteStateSet(initStates, stoichMatrix);
            stateSpace = stateSpace.expand(constraintFunctions, constraintBoundsFinal);
            AfspFull = AfspFull.regenerate(propensities, parameters, stateSpace, constraintCount,speciesNames);
        end

        stateCountOld = stateCount;
        stateCount = stateSpace.getNumStates;
        if useHybrid
            solVec = zeros(stateCount + constraintCount + numODEs, 1);
        else
            solVec = zeros(stateCount + constraintCount, 1);
        end

        solVec(1:stateCountOld) = solVecOld(1:stateCountOld);
        solVec(stateCount+1:end) = solVecOld(stateCountOld+1:end);
    else
        solVec = solVecOld;
    end
end

constraintBoundsFinal = constraintBoundsFinal(1:end-nEscapeSinks);

end

function f = packFspSolution(states, p)
f = ssit.FspVector(states, p);
end

function y = generate_propensity_vector(t, v, x, propensities, parameters)
y = zeros(length(propensities), 1);
wt = propensities{1}.hybridFactorVector(t,parameters,v);
for i = 1:length(propensities)
    if propensities{i}.isFactorizable
        y(i) = wt(i).*propensities{i}.stateDependentFactor(x,parameters);
    else
        y(i) = propensities{i}.hybridJointFactor(t,x,parameters,v);
    end
end
end


