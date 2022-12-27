classdef ExpokitPiecewise
    %EXPOKIT ODE integrator using Expokit.    
    
    properties
        tol = 1.0e-8;
        m = 30;
    end
    
    methods
        function obj = ExpokitPiecewise(m, tol)
        %EXPOKIT Construct an instance of this class
        %   Detailed explanation goes here
        arguments
            m (1,1) double {mustBePositive} = 30
            tol (1,1) double {mustBePositive} = 1.0e-8
        end
        obj.m = m;
        obj.tol = tol;
        end
        
        function [tExport, solutionsNow, fspStopStatus] = solve(obj,...
                tStart, tOut, initSolution, ~, jac, fspErrorCondition)
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
        m = obj.m;
        if ~exist('m','var'); m=15; end
        fspTol = fspErrorCondition.fspTol;
        nSinks = fspErrorCondition.nSinks;
        
        solutionsNow = zeros(size(initSolution,1),length(tOut)-1);
        for iStep = 1:length(tOut)-1
            tStartStep = tOut(iStep);
            tOutStep = tOut(iStep:iStep+1);
            tryAgain=1;
            while tryAgain==1
                [~, ~, ~, tExportStep, solutionsNowStep, ~, tryAgain, te, ye] = ssit.fsp_ode_solvers.mexpv_modified_2(tOutStep(end), jac(tStartStep), initSolution, fspTol/1e5, m,...
                    [], tOutStep, fspTol, [length(initSolution)-nSinks+1:length(initSolution)], tStartStep, fspErrorCondition);
 
                if tryAgain==0;break;end
                if m>300
                    warning('Expokit expansion truncated at 300');
                    [~, ~, ~, tExportStep, solutionsNowStep, ~, tryAgain, te, ye] = ssit.fsp_ode_solvers.mexpv_modified_2(tOutStep(end), jac(tStartStep), initSolution, fspTol/1e5, m,...
                        [], tOutStep, fspTol, [length(initSolution)-nSinks+1:length(initSolution)], tStartStep, fspErrorCondition);
                end
                m=m+5;
            end

            if length(tExportStep)>1&&tExportStep(2)==tOutStep(2)
%                 ikeep = tExportStep(2)==tOutStep(2);
                solutionsNow(:,iStep:iStep+1) = solutionsNowStep';
                initSolution = solutionsNowStep(2, :)';
            else
                iStep=iStep-1;
                break
            end
        end
        if iStep == length(tOut)-1
            tExport = tOut(1:end);
        else
            tExport = tOut(1:iStep+1);
            solutionsNow = solutionsNow(:,1:iStep+1);
        end

        ikeep = ismember(tExport, tOut) & (tExport > tStart);
        solutionsNow = solutionsNow(:,ikeep);
        solutionsNow = mat2cell(solutionsNow, size(solutionsNow,1), ones(1, size(solutionsNow,2)));
        
        if (~isempty(te))
            sinks = ye(end+1-fspErrorCondition.nSinks:end);
            te = max(te,fspErrorCondition.tInit+(tOut(2)-fspErrorCondition.tInit)/1e8);
            errorBound = fspErrorCondition.fspTol*(te-fspErrorCondition.tInit)/(fspErrorCondition.tFinal-fspErrorCondition.tInit);
            if te==max(tOut)
                err = sum(sinks);
            else
                err = sum(sinks*fspErrorCondition.nSinks);
            end
            if err>=errorBound
                fspStopStatus = struct('i_expand', true, ...
                    't_break', te, ...
                    'sinks', sinks, ...
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

