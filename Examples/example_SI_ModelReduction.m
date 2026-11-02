%% SSIT/Examples/example_SI_ModelReduction
%% Preliminaries:
% clear
% close all
% addpath(genpath('../src'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Reduction
%   * Create reduced FSP models using different types of transformations
%     and compare their accuracy and speed to the full FSP solution.
%   * Both time-invariant and time-varying models are supported. For
%     time-varying models, each time-dependent propensity should factorize
%     as f(t)*g(x) (see model 4 below); otherwise the full FSP matrix must
%     be rebuilt and projected at every time step, which removes any
%     speedup from the reduction.
%   * Model reduction is not available for hybrid (ODE/FSP) models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First, choose a model on which to illustrate the reduction approximation,
% or you can create your own. Here are the example options defined below:
%       (1)   2-species Poisson process
%       (2)   Time-varying 2-species Poisson process
%       (3)   3-species repressilator model
%       (4)   Time-varying 3-species repressilator model
%       (5)   Time-varying 4-state STL1 model (loaded from ExampleSaveFiles)
%       (6)   Time-varying bursting transcription/translation (TXTL) model
%       (7)   Time-varying 2-species toggle switch

% Every combination of the models in 'testModels' and the reductions in
% 'reductionTypes' will be tested, and the results summarized in a table.
% 'modelNames' and 'isTimeVarying' are indexed by model number, and
% 'allowsTimeVarying' is indexed by position in 'reductionTypes'.
testModels = [1:7];
modelNames = {'2-Species Poisson','Time-Varying 2-Species Poisson',...
    'Repressilator','TV Repressilator',...
    '4-State STL1','Time-Varying Bursting TXTL','Time-Varying Toggle Switch'};
isTimeVarying = [false,true,false,true,true,true,true];

reductionTypes = {'POD2',...
    'Log Lump QSSA',...
    'Linear State Lumping',...
    'Logarithmic State Lumping',...
    'QSSA',...
    'Eigen Decomposition'};  
allowsTimeVarying = [true,false,true,true,false,false];

reductionOrder = 30;
podTimeSetSize = 200;     % only needed for the POD reduction scheme
makePlots = false;        % for each model, plot the full FSP (thick black) and all reductions (colors) together
% Combinations of time-varying models (isTimeVarying) with reductions that
% are not appropriate for time-varying problems (allowsTimeVarying = false)
% are skipped.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose which type of model reduction to apply.
% Options that can be used for time-invariant AND time-varying models:
%   'Linear State Lumping' (LNSL) --
%       State space is divided into linearly distributed bins.
%       Then within each bin, the proability distribution is assumed to
%       be constant. The number of bins for each species is defined by
%       reductionOrder.
%   'Logarithmic State Lumping' (LGSL) --
%       State space is divided into logarithmically distributed bins.
%       Then within each bin, the proability distribution is assumed to
%       be constant. The number of bins for each species is defined by
%       reductionOrder.
%   'Proper Orthogonal Decomposition' (POD) -
%       The infinitesimal generator is projected onto the orthonormal
%       basis spanned by a previous solution of the full master
%       equation at discrete time points.  This is done by performing
%       SVD on the solutions and choosing the output space
%       corresponding to the 'reductionOrder' largest singular values.
%       For best use, this reduction should be found using a
%       fine time resolution in the calculation of the FSP (set by
%       'podTimeSetSize'). Because the POD requires a solution of the
%       FSP, this reduction is usually only helpful for situations where
%       many solutions are needed (e.g., during model fitting).
%   'POD2' (or 'POD 2nd') --
%       Extension of POD that also includes the time derivatives of the CME
%       solution (i.e., A*P(t_i)) in the set of vectors onto which the CME
%       is projected. Only the vectors corresponding to the top
%       'reductionOrder' singular values are kept for the projection.
%   'No Transform' -
%       Test case where no reduction is applied.
%
% Options that are only applicable to time-invariant models (set
% 'allowsTimeVarying' to false for these). Their bases are built from the
% generator at the initial time, so they do not track time-varying rates:
%   'QSSA' --
%       Species listed in 'qssaSpecies' are assumed to be at
%       quasi-steady state given the other species.
%   'Log Lump QSSA' (LGQSSA) --
%       Same as LGSL, but where the distribution within each lump is
%       assumed to be distributed according to the quasi-steady state
%       assumption.
%   'Eigen Decomposition' (ED) --
%       The infinitesimal generator is projected onto the eigenvectors
%       corresponding to its least negative eigenvalues. The number of
%       eigenvectors for this projection is defined by 'reductionOrder'.
%
% See src/+ssit/+fsp_model_reduction/getTransformMatrices.m for
% additional, less-tested options (e.g., 'Dynamic Mode Decomposition',
% 'POD Update').
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Results are saved after each combination to a separate table for each
% model (named 'results_<modelName>'), using '-append' so that the results
% for one model can be updated without re-running the others.
resultsFile = 'example_SI_ModelReductionResults.mat';
% One color per reduction type, so each method has the same color in
% every model's plots.
redColors = lines(length(reductionTypes));

for iModel = 1:length(testModels)
    testModel = testModels(iModel);

    % Define SSIT Model
    % SSIT models are defined as usual:
    switch testModel
        case 1 % 2-Species Poisson
            Model1 = SSIT;
            Model1.species = {'x1';'x2'};
            Model1.initialCondition = [0;0];
            Model1.propensityFunctions = {'kr';'gr*x1';'kp';'gp*x2'};
            Model1.stoichiometry = [1,-1,0,0;0,0,1,-1];
            Model1.parameters = ({'kr',40;'gr',1;'kp',20;'gp',1});
            Model1.fspOptions.initApproxSS = false;
            Model1.tSpan = linspace(0,5,11);
            indTimes = [1,4,7,11];
            indSpecies = [2];
            qssaSpecies = [1];  % Species that are replaced with QSSA
        case 2 % Time-Varying 2-Species Poisson
            Model1 = SSIT;
            Model1.species = {'x1';'x2'};
            Model1.initialCondition = [0;0];
            Model1.propensityFunctions = {'kr';'kr*cos(omega*t)';'gr*x1';'kp';'gp*x2'};
            Model1.stoichiometry = [1,1,-1,0,0;0,0,0,1,-1];
            Model1.parameters = ({'kr',40;'gr',1;'kp',20;'gp',1;'omega',2*pi});
            Model1.fspOptions.initApproxSS = false;
            Model1.tSpan = linspace(0,10,11);
            indTimes = [1,4,7,11];
            indSpecies = [2];
            qssaSpecies = [1];  % Species that are replaced with QSSA
        case 3 % 3-Species Repressilator Model
            Model1 = SSIT('Repressilator');
            Model1.tSpan = linspace(0,20,11);
            indTimes = [1,4,7,11];
            indSpecies = [3];
            qssaSpecies = [1,2];  % Species that are replaced with QSSA
        case 4 % Time-varying 3-Species Repressilator Model
            Model1 = SSIT('Repressilator');
            % Split the time-varying production of x1 into two reactions so
            % that each propensity factorizes as f(t)*g(x). A non-factorizable
            % propensity (e.g., f(t)+g(x)) forces the full FSP matrix to be
            % rebuilt and projected at every time step, which removes any
            % speedup from the reduced model.
            Model1.propensityFunctions{1} = 'kn1*(1/(1+a*(x2^n)))';
            Model1.propensityFunctions{end+1} = 'kn0*(1+cos(2*pi*t))';
            Model1.stoichiometry(:,end+1) = [1;0;0];
            Model1.tSpan = linspace(0,20,11);
            indTimes = [1,4,7,11];
            indSpecies = [3];
            qssaSpecies = [1,2];  % Species that are replaced with QSSA
        case 5 % Time-varying 4-state STL1 Model
            % Load pre-computed FSP solutions:
            S = load('ExampleSaveFiles/example_11_LoadingandFittingData_MH.mat','STL1_4state');
            % Make a copy of the STL1 model to set up for model reduction:
            Model1 = S.STL1_4state;       
            % 4-state STL1:
            % Set up to solve FSP solution again to time following expansion:
            Model1.fspOptions.initApproxSS = false;
            Model1.tSpan = linspace(0,50,11);
            indTimes = [1,4,7,11];
            indSpecies = [5];
            qssaSpecies = [1,2,3,4];  % Species that are replaced with QSSA
        case 6 % Time-varying Bursting TXTL Model
            % Gene switches on/off with a periodically modulated on-rate;
            % active genes transcribe RNA, which is translated to protein.
            Model1 = SSIT();
            Model1.species = {'goff','gon','rna','prot'};
            Model1.initialCondition = [2;0;0;0];
            Model1.propensityFunctions = {'kon*goff*(1+sin(2*pi*t))';'koff*gon';'kr*gon';'gr*rna';'kp*rna';'gp*prot'};
            Model1.stoichiometry = [-1,1,0,0,0,0;1,-1,0,0,0,0;0,0,1,-1,0,0;0,0,0,0,1,-1];
            Model1.parameters = ({'kon',0.5';'koff',1;'kr',20;'gr',1;'kp',5;'gp',1});
            Model1.tSpan = linspace(0,5,201);
            indTimes = [51,101,151,201];
            indSpecies = [4];
            qssaSpecies = [3];  % Species that are replaced with QSSA
        case 7 % Time-Varying Toggle Switch
            % Mutually repressing LacI/LamCI with an extra LamCI
            % degradation reaction.
            Model1 = SSIT();
            Model1.parameters = {'kb',10;'ka',80;'M',20;'g',1};
            Model1.species = {'LacI';'LamCI'};
            Model1.stoichiometry = [1,-1,0, 0, 0;...
                                    0, 0,1,-1,-1];
            Model1.propensityFunctions = {'kb+ka*M^3/(M^3+LamCI^3)';...
                'g*LacI';...
                'kb+ka*M^3/(M^3+LacI^3)';...
                'g*LamCI';...
                'g*(1+cos(2*pi/10*t))*LamCI';...
                };
            Model1.initialCondition = [0;0];
            Model1.customConstraintFuns = {'(LacI-3).^2.*(LamCI-3).^2'};
            Model1.tSpan = linspace(0,100,11);
            indTimes = [1,4,7,11];
            indSpecies = [2];
            qssaSpecies = [1];  % Species that are replaced with QSSA

    end

    % Solve the original Model (once per model, for comparison)
    % Solve once to get the necessary FSP projection space.
    Model1 = Model1.formPropensitiesGeneral(['ModRedModel',num2str(testModel)]);
    % Check that the time-varying flag matches the model's propensities.
    if any(cellfun(@(f)f.isTimeDependent,Model1.propensitiesGeneral))~=isTimeVarying(testModel)
        warning('isTimeVarying(%d) does not match the propensities of model ''%s''.',...
            testModel,modelNames{testModel})
    end
    Model1 = Model1.solve(solver='FSP');

    % Solve again to record FSP solution time following expansion.
    tic
    Model1 = Model1.solve(solver='FSP');
    fullModelSolveTime = toc;
    fspSoln = Model1.Solutions;

    % Turn off further FSP expansion.
    Model1.fspOptions.fspTol = inf;

    % If using the POD, we will also need to generate basis set using solution
    % at finer resolution. Note -- this means that the POD will be inefficient
    % for the initial set up of the reduction.  The benefits typically come
    % from solving the model multiple times with different parameters sets.
    % This finer solution is computed only once per model and is shared by
    % all POD-type reductions.
    if any(ismember(reductionTypes,{'Proper Orthogonal Decomposition','POD','POD 2nd','POD2'}))
        tSpan = Model1.tSpan;
        Model1.tSpan = linspace(min(Model1.tSpan),max(Model1.tSpan),podTimeSetSize);
        fspSolnFine = Model1.solve(solver='fsp',returnType='soln');
        Model1.tSpan = tSpan;
    end

    nSpecies = length(Model1.species);

    % Set up comparison figures for this model: one figure of means +/- SD
    % versus time, and one figure of marginal distributions at indTimes for
    % each species in indSpecies. Figures for model iModel are numbered
    % 100*iModel+1 (means) and 100*iModel+1+j (marginals of species j).
    % The full FSP solution is plotted in thick black, and each reduction
    % is added below in its own color. On the means plot, each curve is
    % shifted slightly in time (by xShift, about 1% of the time range per
    % curve, centered on the true times) so that the error bars do not
    % overlap.
    if makePlots
        xStep = 0.01*(max(Model1.tSpan)-min(Model1.tSpan));
        xShift = @(k)xStep*(k-length(reductionTypes)/2);  % k = 0 (full), iRed (reductions)
        fMeans = figure(100*iModel+1); clf(fMeans);
        sgtitle(fMeans,[modelNames{testModel},' - means \pm SD']);
        fMarg = gobjects(1,length(indSpecies));
        for j = 1:length(indSpecies)
            fMarg(j) = figure(100*iModel+1+j); clf(fMarg(j));
            sgtitle(fMarg(j),[modelNames{testModel},' - ',...
                Model1.species{indSpecies(j)},' marginal distributions']);
        end
        plotReductionComparison(fMeans,fMarg,fspSoln,Model1,indSpecies,indTimes,...
            'Full FSP',{'Color','k','LineWidth',3},xShift(0));
    end

    % Load any previously saved results for this model.
    resultsVarName = matlab.lang.makeValidName(['results_',modelNames{testModel}]);
    if isfile(resultsFile)&&ismember(resultsVarName,who('-file',resultsFile))
        tmp = load(resultsFile,resultsVarName);
        modelResults = tmp.(resultsVarName);
    else
        modelResults = table();
    end

    for iRed = 1:length(reductionTypes)
        reductionType = reductionTypes{iRed};
        status = 'OK';
        redModelSolveTime = NaN;
        maxKS = NaN;

        if isTimeVarying(testModel)&&~allowsTimeVarying(iRed)
            status = 'Skipped: reduction requires a time-invariant model';
        elseif strcmp(reductionType,'QSSA')&&isempty(qssaSpecies)
            status = 'QSSA Reduction Requires at least one QSSA Species';
        elseif strcmp(reductionType,'QSSA')&&length(Model1.species)==length(qssaSpecies)
            status = 'QSSA Reduction Requires at least one non-QSSA Species';
        else
            if ismember(reductionType,{'Proper Orthogonal Decomposition','POD','POD 2nd','POD2'})
                fspSoln2 = fspSolnFine;
            else
                fspSoln2 = fspSoln;
            end

            % try
                % Solving the reduced models
                % Make a copy of the full model.
                Model2 = Model1;

                % Set the solver to use ModelReduction. FSP expansion is
                % already suppressed (fspTol = inf above), as required when
                % using model reductions.
                Model2.modelReductionOptions.useModReduction = true;

                % Set type and order of Model Reduction
                Model2.modelReductionOptions.reductionType = reductionType;
                Model2.modelReductionOptions.reductionOrder = reductionOrder;
                Model2.modelReductionOptions.qssaSpecies = qssaSpecies;

                % Call SSIT to compute the model reduction transformation matrices.
                % TODO -- This has been seen to fail on "Model 4" with Linear
                % State Lumping. The model numbering has since changed, so
                % recheck which model (STL1 is now model 5) is affected.
                Model2 = Model2.computeModelReductionTransformMatrices(fspSoln2);

                % Solve the reduced model once to construct the reduced FSP
                % matrices, then solve again to record the reduced solution time.
                Model2 = Model2.solve(solver='fsp');
                tic
                Model2 = Model2.solve(solver='fsp');
                redModelSolveTime = toc;
                fspSolnRed = Model2.Solutions;

                % Compare results
                % Kolmogorov-Smirnov distance between the full and reduced
                % marginal distributions at the final time.
                finalError = zeros(1,nSpecies);
                for i = 1:nSpecies
                    finalError(i) = max(abs(squeeze(cumsum(sum(double(fspSoln.fsp{end}.p.data),setdiff((1:nSpecies),i)))) - ...
                        squeeze(cumsum(sum(double(fspSolnRed.fsp{end}.p.data),setdiff((1:nSpecies),i))))));
                end
                maxKS = max(finalError);
            % catch err
            %     status = err.message;
            % end
        end

        % Update this model's results table (replacing any earlier result
        % for the same reduction type) and save it to the results file.
        newRow = table({modelNames{testModel}},{reductionType},reductionOrder,...
            fullModelSolveTime,redModelSolveTime,...
            fullModelSolveTime/redModelSolveTime,maxKS,{status},...
            'VariableNames',{'Model','Reduction','ReductionOrder',...
            'FullTime','ReducedTime','SpeedUp','MaxKS','Status'});
        if ~isempty(modelResults)
            modelResults(strcmp(modelResults.Reduction,reductionType),:) = [];
        end
        modelResults = [modelResults;newRow];
        saveStruct = struct(resultsVarName,modelResults);
        if isfile(resultsFile)
            save(resultsFile,'-struct','saveStruct','-append');
        else
            save(resultsFile,'-struct','saveStruct');
        end
        % disp(newRow)

        if makePlots&&strcmp(status,'OK')
            % Add this reduction to the model's comparison figures as a
            % dashed colored line.
            plotReductionComparison(fMeans,fMarg,fspSolnRed,Model1,indSpecies,indTimes,...
                reductionType,{'Color',redColors(iRed,:),'LineWidth',1.5,'LineStyle','--'},xShift(iRed));
        end
        %% Summary table of errors and efficiencies
        % Combine the saved results for all models (including models run in
        % earlier sessions).        
        resultsVars = who('-file',resultsFile,'results_*');
        savedResults = struct2cell(load(resultsFile,resultsVars{:}));
        resultsTable = vertcat(savedResults{:});
        disp(resultsTable)
    end
end

%% Local functions
function plotReductionComparison(fMeans,fMarg,soln,Model,indSpecies,indTimes,label,lineProps,xShift)
% Add one FSP solution (full or reduced) to a model's comparison figures:
%   fMeans -- mean (line) with +/- SD error bars versus time, with one
%             subplot per species in indSpecies (square grid). The curve is
%             shifted by xShift in time so that error bars from different
%             reductions do not overlap. Error bars are drawn at no more
%             than ~20 time points to keep long time series readable.
%   fMarg  -- array of figures, one per species in indSpecies, each
%             showing the marginal distributions at the times in indTimes
%             in a square grid (e.g., 2x2 for 4 times, 3x3 for 9 times).
% 'label' is used in the legend, and 'lineProps' sets color/width/style.
nSpecies = length(Model.species);
nTimes = length(soln.fsp);
nPlotSp = length(indSpecies);
nRows = ceil(sqrt(nPlotSp));
nPlotT = length(indTimes);
nRowsT = ceil(sqrt(nPlotT));
t = Model.tSpan(1:nTimes)+xShift;
indErr = unique(round(linspace(1,nTimes,min(nTimes,20))));
lineColor = lineProps{find(strcmp(lineProps,'Color'),1)+1};
for j = 1:nPlotSp
    s = indSpecies(j);

    % Means and standard deviations versus time.
    mu = zeros(1,nTimes);
    sd = zeros(1,nTimes);
    for it = 1:nTimes
        pmf = getMarginal(soln.fsp{it},s,nSpecies);
        x = 0:length(pmf)-1;
        mu(it) = x*pmf;
        sd(it) = sqrt(max(0,(x.^2)*pmf-mu(it)^2));
    end
    set(0,'CurrentFigure',fMeans);
    ax = subplot(nRows,nRows,j); hold(ax,'on');
    plot(ax,t,mu,lineProps{:},'DisplayName',label);
    errorbar(ax,t(indErr),mu(indErr),sd(indErr),'LineStyle','none',...
        'Color',lineColor,'LineWidth',1,'HandleVisibility','off');
    xlabel(ax,'Time'); ylabel(ax,[Model.species{s},' count']);
    if j==1, legend(ax,'Location','best'); end

    % Marginal distributions at the selected times.
    set(0,'CurrentFigure',fMarg(j));
    for k = 1:nPlotT
        pmf = getMarginal(soln.fsp{indTimes(k)},s,nSpecies);
        ax = subplot(nRowsT,nRowsT,k); hold(ax,'on');
        stairs(ax,(0:length(pmf)-1)-0.5,pmf,lineProps{:},'DisplayName',label);
        title(ax,sprintf('t = %.3g',Model.tSpan(indTimes(k))));
        xlabel(ax,[Model.species{s},' count']); ylabel(ax,'Probability');
        if k==1, legend(ax,'Location','best'); end
    end
end
end

function pmf = getMarginal(fspTimePoint,s,nSpecies)
% Marginal distribution (column vector) of species s at one FSP time point.
pmf = sum(double(fspTimePoint.p.data),setdiff(1:nSpecies,s));
pmf = pmf(:);
end
