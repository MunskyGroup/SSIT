function app = runFsp(app,fsp_soln_style)
% This function runs the FSP calculation after the "Run FSP" Button is
% selected within the GUI.
if nargin<2
    fsp_soln_style = 'default';
end
%% Provides an indicator when the FSP analysis is running for the user
if verLessThan('matlab','9.4') % This if statement allows for the script to be run in other versions than the latest version
    app.FspRunButton.BackgroundColor = 'yellow';
    pause(2)
else
    try
%         d_prog_bar = waitbar(0,'Running FSP Computation');
        d_prog_bar = [];%waitbar(0,'Running FSP Computation');
    catch
        d_prog_bar=[];
    end
end
%%
switch app.FspIntegratorTolerance.Value
    case 'Relaxed'
        relTol = 1e-2; absTol = 1e-4;
    case 'Moderate'
        relTol = 1e-3; absTol = 1e-6;
    case 'Strict'
        relTol = 1e-4; absTol = 1e-12;
end
%%
x0 = eval(app.FspInitCondField.Value)';   % Pulls the inital conditions specified
stoichMatrix = app.ReactionsTabOutputs.stoichMatrix;
if (size(stoichMatrix, 1) ~= size(x0, 1))
    stoichMatrix = stoichMatrix';
    app.ReactionsTabOutputs.stoichMatrix=stoichMatrix;
end

parsDict = containers.Map(app.ReactionsTabOutputs.parameters(:,1), app.ReactionsTabOutputs.parameters(:,2));
inputs = app.ReactionsTabOutputs.inputs;

propstrings = ssit.SrnModel.processPropensityStrings(app.ReactionsTabOutputs.propensities,inputs,parsDict,"str",app.ReactionsTabOutputs.varNames);

% Note - this needs to be cleaned up.  I would like to separate expressions
% so that 'x' terms are factorized out if possible.
% propstrings{1} = str2sym('k');
% propstrings{2} = str2sym('g*x1');
% propensities{i} = ssit.Propensity.createFromSym(propstrings{i}, stoichMatrix(:,i), i);

n_reactions = length(app.ReactionsTabOutputs.propensities);
% propensities = cell(n_reactions, 1);
% for i = 1:n_reactions
%     propstrings{i} = strrep(propstrings{i},'..','.');
%     propensities{i} = ssit.Propensity.createFromString(propstrings{i}, stoichMatrix(:,i), i);
% end


if isfield(app.FspTabOutputs,'PropensitiesGeneral')
    PropensitiesGeneral = app.FspTabOutputs.PropensitiesGeneral;
else

    sm = cell(1,n_reactions);
    logicTerms = cell(1,n_reactions);
    logCounter = 0;
    for i = 1:n_reactions
        st = app.ReactionsTabOutputs.propensities{i};
        for jI = 1:size(app.ReactionsTabOutputs.inputs,1)
            st = strrep(st,app.ReactionsTabOutputs.inputs{jI,1},['(',app.ReactionsTabOutputs.inputs{jI,2},')']);
        end
        [st,logicTerms{i},logCounter] = ssit.Propensity.stripLogicals(st,app.ReactionsTabOutputs.varNames,logCounter);
        sm{i} = str2sym(st);
    end

    % if obj.useHybrid
    %     PropensitiesGeneral = ssit.Propensity.createAsHybridVec(sm, obj.stoichiometry,...
    %         obj.parameters, obj.species, obj.hybridOptions.upstreamODEs, logicTerms, prefixName);
    % else
    PropensitiesGeneral = ssit.Propensity.createAsHybridVec(sm, ...
        app.ReactionsTabOutputs.stoichMatrix,...
        app.ReactionsTabOutputs.parameters, ...
        app.ReactionsTabOutputs.varNames, ...
        [], logicTerms, 'TEMPPropens');
    app.FspTabOutputs.PropensitiesGeneral = PropensitiesGeneral;
end
% end


% app.FspTabOutputs.propensities = propensities;

if (size(x0, 2) > 1)
    x0 = x0';
end

useMex = 0;%app.FspUseMexCheckBox.Value;
tSpan = unique(eval(app.FspPrintTimesField.Value));        % Pulls the time array from the app
fspTol = app.FspErrorTolField.Value;            % Pulls the specified error tolerance for the FSP Analysis
fConstraints = readConstraintsForAdaptiveFsp(app);
bConstraints = app.FspTabOutputs.bounds';
%%
isTimeInvar = size(app.ModelInputTable.Data,1)==0;
if isTimeInvar
    useExpokit = 1;
else
    useExpokit = 0;
end
   
app.FspRunningStatus.Text = 'FSP is running...';
tic
verbose = 1;
odeSolver = 'auto';

[solutions, bConstraints, app.FspTabOutputs.stateSpace] = ssit.fsp.adaptiveFspSolve(  tSpan,...
                                                        x0,...
                                                        [1.0], ...
                                                        stoichMatrix, ...
                                                        PropensitiesGeneral, ...  % need change and add parameters
                                                        [app.ReactionsTabOutputs.parameters{:,2}]',...
                                                        fspTol, ...
                                                        fConstraints, ...
                                                        bConstraints,...
                                                        verbose, ...
                                                        relTol, ...
                                                        absTol, ...
                                                        odeSolver, ...
                                                        app.FspTabOutputs.stateSpace,...
                                                        app.FspPiecewiseCheckBox.Value,...
                                                        app.initApproxSS.Value,...
                                                        app.ReactionsTabOutputs.varNames);

%% Re-sets the indicator for running of FSP
app.FspTabOutputs.bounds = bConstraints';
app.FspTabOutputs.solutions = solutions;
app.FspRunningStatus.Text = ['FSP completed in ',num2str(toc,3),'s'];
if verLessThan('matlab','9.4') % THis if statement allows for user to know it is ready to go for any version of Matlab
    app.FspRunButton.BackgroundColor = 'green';
    pause(2)
else
    delete(d_prog_bar); % Closes the Running Dialog Box
end
app.FspConstraintTable.Data(:,3) = num2cell(app.FspTabOutputs.bounds');   % Sets the bounds from constraints to the ones found in the FSP analysis
end
