classdef OdeSuite < ssit.fsp_ode_solvers.OdeSolver
    %ODESUITE ODE integrator using MATLAB's ODE suite.        
    properties
        relTol (1,1) double {mustBePositive} = 1.0e-4
        absTol (1,1) double {mustBePositive} = 1.0e-8   
        solver = @ode23s
        maxStep (1,1) double {mustBePositive} = 0.01
%         numODEs = 0
    end
    
    methods
        function obj = OdeSuite(relTol, absTol, solver, maxStep)
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
            solver = 'ode23s'
            maxStep (1,1) double {mustBePositive} = 1.0
%             numODEs = 0
        end
        
        obj.relTol = relTol;
        obj.absTol = absTol;
        obj.solver = str2func(solver);
        obj.maxStep = maxStep;
%         obj.numODEs = numODEs;
        end
        
        function [tExport, solutionsNow, fspStopStatus] = solve(obj, tStart,...
                tOut, initSolution, rhs, jac, fspErrorCondition, fixedEvents)
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
        odeEvent = @(t,p) fspOdesuiteEvent(t, p, fspErrorCondition);  
        if length(tOut)>1
            maxStep = min(obj.maxStep,min(tOut(2:end)-tOut(1:end-1))/2);
        else
            maxStep = min(obj.maxStep,(tOut-tStart)/2);
        end
        if ~isempty(jac)
            ode_opts = odeset('Events', odeEvent, 'Jacobian', jac, 'relTol',...
                obj.relTol, 'absTol', obj.absTol,'Vectorized','on',...
                'MaxStep',maxStep,'JPattern',jac(rand,rand(size(initSolution))~=0));
        else
            ode_opts = odeset('Events', odeEvent, 'relTol',obj.relTol,...
                'absTol', obj.absTol,'Vectorized','off','MaxStep',maxStep);
        end
               
        % Initialize output storage
        tExportAll = [];
        solutionsAll = [];
        tCurrentStart = tStart;
        currentSolution = initSolution;
        
        isStub = true;
        while isStub

            % Identify fixed event times within integration interval
            fixedEventTimesInInterval = [];
            if ~isempty(fixedEvents) && isfield(fixedEvents, 'times')
                fixedEventTimesInInterval = fixedEvents.times(fixedEvents.times > tCurrentStart & fixedEvents.times < tOut(end));
            end

            % Integrate in segments between fixed events
            if ~isempty(fixedEventTimesInInterval)&&fixedEventTimesInInterval(1)<tOut(end)
                tSpan = sort(unique([tCurrentStart; tOut; fixedEventTimesInInterval']));
                tSpan = tSpan(tSpan<=fixedEventTimesInInterval(1)&tSpan<=tOut(end));
            else
                % No fixed events: integrate normally
                tSpan = sort(unique([tCurrentStart; tOut]));
            end
            tSpan = tSpan(tSpan>=tCurrentStart);

            if length(tSpan)==2
                tSpan = [tSpan(1),mean(tSpan),tSpan(2)];
            end
            [tExport, solutionsNow, te, ye, ~] =  obj.solver(rhs, tSpan, currentSolution, ode_opts);

            if ~isempty(te)&&(length(tExport)<2||te(end)<tSpan(2))
                solutionsAll = {};
                tExportAll = tExport(1);
                break
            else

                % Apply fixed time special event if applicable
                if ~isempty(fixedEventTimesInInterval)
                    if tExport(end) == fixedEventTimesInInterval(1)
                        fixedEventIdx = find(abs(fixedEventTimesInInterval(1) - fixedEvents.times) < 1e-10, 1, 'first');
                        if ~isempty(fixedEventIdx)
                            matrixIdx = fixedEvents.matrixInds(fixedEventIdx);
                            currentSolution = fixedEvents.matrices{matrixIdx} * solutionsNow(end, :)';
                            solutionsNow(end, :) = currentSolution';
                        end
                    end
                end
                tExportAll = [tExportAll;tExport];
                solutionsAll = [solutionsAll;solutionsNow];
            end

            if isempty(fixedEventTimesInInterval)||(tExport(end)>=tOut(2))
                isStub = false;
            else
                isStub = true;
                % Start next segment at slightly later time (just after
                % fixed event)
                tCurrentStart = tExport(end)*(1+1e-10);
            end
        end

        % Sort to return results only at requested times.
        ikeep = ismember(tExportAll, tOut) & (tExportAll > tStart);
        if ~isempty(ikeep)&&sum(ikeep)>0
            tExportAll = tExportAll(ikeep);
            solutionsAll = solutionsAll(ikeep,:)';
        end

        tExport = tExportAll;
        solutionsNow = solutionsAll;
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
% error_bound = fspErrorCheck.fspTol*(t-fspErrorCheck.tInit)/(fspErrorCheck.tFinal-fspErrorCheck.tInit);

if isinf(fspErrorCheck.fspTol)
    error_bound = inf;
else
    error_bound = fspErrorCheck.fspTol*...
        (t-fspErrorCheck.tInit)/(fspErrorCheck.tFinal-fspErrorCheck.tInit);
end

% val = max(sinks)*(fspErrorCheck.nSinks-fspErrorCheck.nEscapeSinks) - error_bound;
val = sum(sinks) - error_bound;
% val is used in ode23s.  It is an indicator if the
% solution should continue or be interupted. 
% When val crosses zero to become positive (i.e. if the
% error becomes too large), this will trigger the solution scheme to execute the
% function
% See example here:
% https://www.mathworks.com/help/matlab/ref/odeevent.html

terminal = true; % if the event is triggered, then it will be terminal.
direction = 1; % check to see if the error is increasing.
end
