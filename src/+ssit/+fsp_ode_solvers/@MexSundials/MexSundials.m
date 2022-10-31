classdef MexSundials < ssit.fsp_ode_solvers.OdeSolver
    %MEXSUNDIALS ODE integrator using MEX interface to SUNDIALS.    
    
    properties
        relTol (1,1) double {mustBePositive} = 1.0e-4
        absTol (1,1) double {mustBePositive} = 1.0e-8        
    end
    
    methods
        function obj = MexSundials(relTol, absTol)
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
        arguments
            relTol (1,1) double {mustBePositive} = 1.0e-4
            absTol (1,1) double {mustBePositive} = 1.0e-8
        end
        
        obj.relTol = relTol;
        obj.absTol = absTol;
        end
        
        function [tExport, solutionsNow, fspStopStatus] = solve(obj, tStart, tOut, initSolution, rhs, jac, fspErrorCondition)
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
        fOut = @(t,p) p;        
        [solutionsNow, fspStopStatus] = ...
            FspCVodeMex(tStart, tOut,...
            rhs, ...
            initSolution,...
            fOut,...
            @fspErrorChecking,...
            fspErrorCondition);
        tExport = tOut(1:size(solutionsNow, 2));
        end
    end
end

function [stopHere, fspErrorStatus] = fspErrorChecking(t, p, fspErrorCheck)
sinks = p(end+1-fspErrorCheck.nSinks:end);
error_bound = fspErrorCheck.fspTol*t/fspErrorCheck.tFinal;
if (max(sinks)*fspErrorCheck.nSinks > error_bound)
    stopHere = double(1);
    fspErrorStatus = struct('i_expand', true, 't_break', t, 'sinks', sinks, 'error_bound', error_bound);
else
    stopHere = double(0);
    fspErrorStatus = struct('i_expand', false, 't_break', t, 'sinks', sinks, 'error_bound', error_bound);
end
end
