classdef Expokit
    %EXPOKIT ODE integrator using Expokit.    
    
    properties
        tol = 1.0e-8;
        m = 30;
        version = 'mexpv_modified_2';
    end
    
    methods
        function obj = Expokit(m, tol, version)
        %EXPOKIT Construct an instance of this class
        %   Detailed explanation goes here
        arguments
            m (1,1) double {mustBePositive} = 30
            tol (1,1) double {mustBePositive} = 1.0e-8
            version  = 'mexpv_modified_2'
        end
        obj.m = m;
        obj.tol = tol;
        obj.version = version;
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
        m = min(ceil(size(jac,1)/2),obj.m);
        tryAgain=1;
        if ~exist('m','var'); m=15; end
        fspTol = fspErrorCondition.fspTol;
        nSinks = fspErrorCondition.nSinks;
        
        if strcmp(obj.version,'mexpv_modified_2')
            while tryAgain==1
                SINKS = [length(initSolution)-nSinks+1:length(initSolution)-fspErrorCondition.nEscapeSinks];
                fspSINKS = [length(initSolution)-nSinks+1 : length(initSolution)];
                [~, ~, ~, tExport, solutionsNow, ~, tryAgain, te, ye] = ssit.fsp_ode_solvers.mexpv_modified_2(tOut(end), jac, initSolution, fspTol/1e5, m,...
                    [], tOut, fspTol, SINKS, tStart, fspErrorCondition);
                if tryAgain==0;break;end
                if m>300
                SINKS = [length(initSolution)-nSinks+1:length(initSolution)-fspErrorCondition.nEscapeSinks];
                    warning('Expokit expansion truncated at 300');
                    [~, ~, ~, tExport, solutionsNow, ~, tryAgain, te, ye] = ssit.fsp_ode_solvers.mexpv_modified_2(tOut(end), jac, initSolution, fspTol/1e5, m,...
                        [], tOut, fspTol, SINKS, tStart, fspErrorCondition);
                end
                m=m+5;
            end
        elseif strcmp(obj.version,'expv_modified')
            while tryAgain==1
                SINKS = [];
                [~, ~, ~, tExport, solutionsNow, ~, tryAgain, te, ye] = ...
                    ssit.fsp_ode_solvers.expv_modified(tOut(end), jac, initSolution, 1e-16, m,...
                    [], tOut,fspTol,SINKS,tStart,fspErrorCondition);
                if tryAgain==0;break;end
                if m>300
                    warning('Expokit expansion truncated at 300');
                    [~, ~, ~, tExport, solutionsNow, ~, tryAgain, te, ye] = ...
                        expv_modified(tOut(end), jac, initSolution, fspTol/1e5, m,...
                        [], tOut,fspTol,[],tStart,fspErrorCondition);
                end
                m=m+5;
            end
        end
        
        ikeep = ismember(tExport, tOut) & (tExport > tStart);
        solutionsNow = solutionsNow(ikeep, :)';
        solutionsNow = mat2cell(solutionsNow, size(solutionsNow,1), ones(1, size(solutionsNow,2)));
        
        if (~isempty(te))
            sinks = ye(SINKS);
            errorBound = fspErrorCondition.fspTol*(te-fspErrorCondition.tInit)/(fspErrorCondition.tFinal-fspErrorCondition.tInit);
            if te==max(tOut)
                err = sum(sinks);
            else
                err = sum(sinks*(fspErrorCondition.nSinks-fspErrorCondition.nEscapeSinks));
            end
            if err>=errorBound
                fspStopStatus = struct('i_expand', true, ...
                    't_break', te, ...
                    'sinks', ye(fspSINKS), ...
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

