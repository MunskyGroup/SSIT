function fsp_model_reduction(app,redType,sweepPar,sweepMin,sweepMax)
arguments
    app
    redType
    sweepPar = [];
    sweepMin = [];
    sweepMax = [];
end
phi = [];
%% Find Original Infinite Generator Matrix

% Set up propensity functions
stoichMatrix = app.SSITapp.ReactionsTabOutputs.stoichMatrix;
baseParameters =[app.SSITapp.ReactionsTabOutputs.parameters{:,2}];
parNames = app.SSITapp.ReactionsTabOutputs.parameters(:,1);

if ~isempty(sweepPar)
    nParVals = 10;
    parSets = repmat(baseParameters,nParVals,1);
    sweepParNumber = find(strcmp(parNames,sweepPar));
    parSets(2:nParVals,sweepParNumber) = linspace(sweepMin,sweepMax,nParVals-1);
else
    nParVals=1;
    parSets = baseParameters;
end

for iPar = 1:nParVals
    if isempty(app.fspSoln)||iPar>1
        parsDict = containers.Map(parNames, parSets(iPar,:));
        inputs = app.SSITapp.ReactionsTabOutputs.inputs;

        propstrings = ssit.SrnModel.processPropensityStrings(app.SSITapp.ReactionsTabOutputs.propensities,inputs,parsDict);

        n_reactions = length(app.SSITapp.ReactionsTabOutputs.propensities);
        propensities = cell(n_reactions, 1);
        for i = 1:n_reactions
            propstrings{i} = strrep(propstrings{i},'..','.');
            propensities{i} = ssit.Propensity.createFromString(propstrings{i}, stoichMatrix(:,i), i);
        end

        % Set up the initial state subset
        initStates = eval(app.SSITapp.FspInitCondField.Value)';   % Pulls the inital conditions specified
        stateSpace = ssit.FiniteStateSet(initStates, stoichMatrix);
        fConstraints = readConstraintsForAdaptiveFsp(app.SSITapp);
        stateSpace = stateSpace.expand(fConstraints, app.SSITapp.FspTabOutputs.bounds');
        constraintCount = length(app.SSITapp.FspTabOutputs.fConstraints);
        fspSoln.states = stateSpace.states;

        % Generate the FSP matrix
        fspSoln.Afsp = ssit.FspMatrix(propensities, stateSpace, constraintCount);
        fspSoln.A_total = fspSoln.Afsp.terms{1}.matrix;
        for i = 2:length(fspSoln.Afsp.terms)
            fspSoln.A_total = fspSoln.A_total+fspSoln.Afsp.terms{i}.matrix;
        end
        fspSoln.A_total = fspSoln.A_total(1:end-constraintCount,1:end-constraintCount);

        fspSoln.tOut = unique(eval(app.SSITapp.FspPrintTimesField.Value));        % Pulls the time array from the app

        % Generate initial distribution
        nStates = size(fspSoln.A_total,1);
        fspSoln.P0 = zeros(nStates, 1);
        fspSoln.P0(1:size(initStates,2)) = 1;

        % Solve Original Problem
        tic
        disp('Running full FSP solve.')
        [~, ~, ~, ~, fspSoln.fullSolutionsNow] = ...
            ssit.fsp_ode_solvers.expv_modified(fspSoln.tOut(end), fspSoln.A_total, fspSoln.P0, 1e-8, 30,...
            [], fspSoln.tOut, 1e-3, [], fspSoln.tOut(1));
        tFullSolve = toc;
        app.fspSoln = fspSoln;
        app.fspSoln.tFullSolve = tFullSolve;
    else
        fspSoln = app.fspSoln;
        tFullSolve = app.fspSoln.tFullSolve;
    end

    lastRed = 'Empty';
    for ired=1:length(redType)
        disp(redType{ired})

        switch redType{ired}
            case {'Eigen Decomposition','Eigen Decomposition Initial','Proper Orthogonal Decomposition'...
                    'Balanced Model Truncation (HSV)','Dynamic Mode Decomposition',...
                    'Radial Basis Functions'}
                n = app.NumberofBasisVectorsEditField.Value;
            case {'Linear State Lumping','Logarithmic State Lumping'}
                n = app.GridSizeEditField.Value;
        end

%         try
            tic
            if isempty(phi)||~strcmp(lastRed,redType{ired})
                [phi,phi_inv,redOutputs] = ssit.fsp_model_reduction.getTransformMatrices(redType{ired},n,fspSoln);
            end
            reductionTime(ired) = toc;
            %% Solve the reduced model
            tic
            if strcmp(redType{ired},'Balanced Model Truncation (HSV)')
                sys = ss(fspSoln.A_total,fspSoln.P0,eye(nStates),[]);
                sysred = balred(sys,n,redOutputs.info);
                A_red = sysred.A;
                q0 = sysred.B;
                OutPutC = sysred.C;
            else
                q0 = phi_inv*fspSoln.P0;
                A_red = phi_inv*fspSoln.A_total*phi;
            end
            [~, ~, ~, ~, solutionsNow] = ssit.fsp_ode_solvers.expv_modified(fspSoln.tOut(end), A_red, q0, 1e-8, 30,...
                [], fspSoln.tOut, 1e-3, [], fspSoln.tOut(1));

            if strcmp(redType{ired},'Balanced Model Truncation (HSV)')
                redSolutionsNow = solutionsNow*OutPutC';
            else
                redSolutionsNow = solutionsNow*phi';
                redSolutionsNow = diag(1./sum(redSolutionsNow,2))*redSolutionsNow;
            end
            solutionTime(ired) = toc;

            % Plot the difference over time
            switch app.ErrorTypeDropDown.Value
                case 'One Norm'
                    ErrorVsTime(ired,iPar,:) = sum(abs(fspSoln.fullSolutionsNow-redSolutionsNow),2)';
                case 'KLD'
                    p = max(1e-6,fspSoln.fullSolutionsNow);
                    q = max(1e-6,redSolutionsNow);
                    ErrorVsTime(ired,iPar,:) = 1/2*(sum(p.*log(p./q),2)+sum(q.*log(q./p),2))';
            end

            switch app.plotStyle.Value
                case 'linear'
                    plot(app.errorVsTime,fspSoln.tOut(2:end),squeeze(ErrorVsTime(ired,1,2:end)),'linewidth',2)
                case 'logarithmic'
                    semilogy(app.errorVsTime,fspSoln.tOut(2:end),squeeze(ErrorVsTime(ired,1,2:end)),'linewidth',2)
            end
            lastRed = redType{ired};
%         catch ME
%             disp(['Error Encountered - Skipping to next Method']);
%             ME.message
%             reductionTime(ired) = inf;
%             solutionTime(ired) = inf;
%             ErrorVsTime(ired,iPar,:) = inf*ones(size(fspSoln.tOut));
%         end

    end
end
if ~isempty(sweepPar)
    [~,J] = sort(parSets(:,sweepParNumber));
    plot(app.errorVsPars,parSets(J,sweepParNumber),squeeze(sum(ErrorVsTime(ired,J,:),3))'/length(fspSoln.tOut),'-o','linewidth',3)
else
    app.reductionTimeTable.Data = {};
    app.reductionTimeTable.Data(1,:) = {'Full FSP','N/A',num2str(tFullSolve,3),'1','N/A'};
    for ired=1:length(redType)
        app.reductionTimeTable.Data(ired+1,:) = {redType{ired},num2str(reductionTime(ired),3),...
            num2str(solutionTime(ired),3),num2str(tFullSolve/solutionTime(ired),3),...
            num2str(mean(ErrorVsTime(ired,1,:)),3)};
    end
end