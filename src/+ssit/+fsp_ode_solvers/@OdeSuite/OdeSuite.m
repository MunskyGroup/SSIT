classdef OdeSuite < ssit.fsp_ode_solvers.OdeSolver
    %ODESUITE ODE integrator using MATLAB's ODE suite.        
    properties
        relTol (1,1) double {mustBePositive} = 1.0e-4
        absTol (1,1) double {mustBePositive} = 1.0e-8    
%         numODEs = 0
    end
    
    methods
        function obj = OdeSuite(relTol, absTol)
        % Construct an instance of MexSundials.
        % 
        % Parameters
        % ----------
        %
        % relTol: double
        %   relative tolerance for the solver.
        %
        % absTol: double
        %   absolute tolerance for the solver.
        %
        % numODEs: int
        %   numer of upstream ODEs at the end of the state vector (for
        %   hybrid models.
        arguments
            relTol (1,1) double {mustBePositive} = 1.0e-4
            absTol (1,1) double {mustBePositive} = 1.0e-8
%             numODEs = 0
        end
        
        obj.relTol = relTol;
        obj.absTol = absTol;
%         obj.numODEs = numODEs;
        end
        
        function [tExport, solutionsNow, fspStopStatus] = solve(obj, tStart, tOut,...
                initSolution, rhs, jac, fspErrorCondition, odeSolnScheme)
        %SOLVE Advance the solution of the FSP-truncated CME up until
        %either the final time point or when the FSP error is exceeded for
        %the current set of states.
        %
        % Parameters
        % ----------
        %
        % tStart: double
        %   Start time for the integration.
        %
        % tOut: 1-D array of doubles
        %   Time points to output the solution.
        %
        % initSolution: 1-D array of doubles
        %   Initial solution, including the sink states. These sink states
        %   must be contiguous and found at the end of the vector.
        %
        % rhs: callable
        %   Function to evaluate the action of the time-varying CME matrix.
        %   Callable via the syntax w = rhs(t, v).
        %
        % jac: unused        
        %
        % fspErrorCondition: struct
        %   Information to decide early exit based on FSP error tolerance.
        %   Must have the following fields:
        %   - nSinks: number of sink states.
        %   - fspTol: FSP tolerance.
        %   - tFinal: final time.
        %
        % Returns
        % -------
        %
        % tExport: 1-D array
        %   Time points where solutions are successfully computed.
        %
        % solutionsNow: 2-D array
        %   Exported solutions. Each column corresponds to the solution at
        %   a timepoint. Specifically, ``solutionsNow(:, j)`` is the solution
        %   at time ``tExport(j)``.
        %
        % fspStopStatus: struct
        %   Status when exiting the ODE integrator. Must contain the
        %   following fields:
        %   - needExpansion: (true/false) whether to expand the FSP for
        %   future time steps.
        %   - tExit: exit time.
        %   - sinks: vector of sink probabilities at ``tExit``.
        %   - errorBound: the error bound at ``tExit``. 
        %
        %
        arguments
            obj
            tStart
            tOut
            initSolution
            rhs
            jac
            fspErrorCondition
            odeSolnScheme = 'ode23s';
        end

        odeEvent = @(t,p) fspOdesuiteEvent(t, p, fspErrorCondition);  
        if length(tOut)>1
            maxStep = min(tOut(2:end)-tOut(1:end-1))/2;
        else
            maxStep = (tOut-tStart)/2;
        end
        if ~isempty(jac)
            ode_opts = odeset('Events', odeEvent, 'Jacobian', jac, 'relTol',...
                obj.relTol, 'absTol', obj.absTol,'Vectorized','on',...
                'MaxStep',maxStep);
        else
            ode_opts = odeset('Events', odeEvent, 'relTol',obj.relTol,...
                'absTol', obj.absTol,'Vectorized','off','MaxStep',maxStep);
        end
        tSpan = sort(unique([tStart; tOut]));

        %% Choose ODE solution Scheme
        switch odeSolnScheme
            case 'ode45' % solver based on the fourth and fifth order Runge-Kutta method
                [tExport, solutionsNow, te, ye, ~] =  ode45(rhs, tSpan, initSolution, ode_opts);
            case 'ode23' % similar to ode45 but uses a second and third order Runge-Kutta method
                [tExport, solutionsNow, te, ye, ~] =  ode23(rhs, tSpan, initSolution, ode_opts);
            case 'ode23s' % similar to ode15s but uses a second and third order Runge-Kutta methoddx
                [tExport, solutionsNow, te, ye, ~] =  ode23s(rhs, tSpan, initSolution, ode_opts);
            case 'ode23t' % A trapezoidal rule solver for solving ODEs with a high-index differential-algebraic equation (DAE) form
                [tExport, solutionsNow, te, ye, ~] =  ode23t(rhs, tSpan, initSolution, ode_opts);
            case 'ode23s_SparseWithJPattern'
                sparsityJac = jac(rand)~=0;
                ode_opts = odeset('Events', odeEvent, 'Jacobian', jac,...
                    'JPattern',sparsityJac,...
                    'relTol', obj.relTol, 'absTol', obj.absTol,'Vectorized','on',...
                    'MaxStep',maxStep);
                [tExport, solutionsNow, te, ye, ~] =  ode23s(rhs, tSpan, initSolution, ode_opts);
            case 'ode15s' % Designed for stiff ODEs, it uses a variable-order, variable-step method based on implicit integration formulas
                [tExport, solutionsNow, te, ye, ~] =  ode15s(rhs, tSpan, initSolution, ode_opts);
            case 'ode113' % variable-order solver that uses a combination of Adams-Bashforth and Gear methods
                [tExport, solutionsNow, te, ye, ~] =  ode113(rhs, tSpan, initSolution, ode_opts);
        end 

        if ~isempty(te)&&(length(tExport)<2||te(end)<tExport(2))
            solutionsNow = solutionsNow(1);
            tExport = tExport(1);
        end
        
        ikeep = ismember(tExport, tSpan) & (tExport > tStart);
        solutionsNow = solutionsNow(ikeep, :)';
        solutionsNow = mat2cell(solutionsNow, size(solutionsNow,1), ones(1, size(solutionsNow,2)));
        
        if (~isempty(te))
            SINKS = length(initSolution)-fspErrorCondition.nSinks-fspErrorCondition.numODEs+1:...
                length(initSolution)-fspErrorCondition.nEscapeSinks-fspErrorCondition.numODEs;
            fspSINKS = length(initSolution)-fspErrorCondition.nSinks-fspErrorCondition.numODEs+1 : ...
                length(initSolution)-fspErrorCondition.numODEs;
            sinks = max(0,ye(SINKS));
            errorBound = fspErrorCondition.fspTol*...
                (te-fspErrorCondition.tInit)/(fspErrorCondition.tFinal-fspErrorCondition.tInit);
            if sum(sinks*(fspErrorCondition.nSinks-fspErrorCondition.nEscapeSinks))>=errorBound
                fspStopStatus = struct('i_expand', true, ...
                    't_break', te, ...
                    'sinks', max(0,ye(fspSINKS)), ...
                    'error_bound', errorBound);
            else
                fspStopStatus = struct('i_expand', false,...
                    't_break', [], ...
                    'sinks', [],...
                    'error_bound', []);
            end
        else
            fspStopStatus = struct('i_expand', false,...
                't_break', [], ...
                'sinks', [],...
                'error_bound', []);
        end
        
        end
    end
end


function [val, terminal, direction] = fspOdesuiteEvent(t, p, fspErrorCheck)
arguments
    t
    p
    fspErrorCheck
end
sinks = p(end-fspErrorCheck.nSinks-fspErrorCheck.numODEs+1:end-fspErrorCheck.numODEs-fspErrorCheck.nEscapeSinks);
error_bound = fspErrorCheck.fspTol*(t-fspErrorCheck.tInit)/(fspErrorCheck.tFinal-fspErrorCheck.tInit);

val = max(sinks)*(fspErrorCheck.nSinks-fspErrorCheck.nEscapeSinks) - error_bound;
terminal = 1;
direction = 1;
end
