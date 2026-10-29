%% SSIT/Examples/example_SI_ModelReduction
%% Preliminaries:
% clear
% close all
% addpath(genpath('../src'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Reduction
%   * Create reduced FSP models using different types of transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First, choose a model on which to illustrate the reduction approximation,
% or you can create your own. Here are the example options defined below:
%       (1)   2-species Poisson process
%       (2)   Constant 3-Species repressilator model
%       (3)   Time-varying 3-species repressilator model
%       (4)   Time-varying 4-state STL1 Model
%       (5)   Time-varying TXTL model from tests/modelReductionTest.m
%       (6)   Goutsias gene regulation model (from example_Benchmarking.m)
%       (7)   2-species toggle switch with gene states (from example_Benchmarking.m)
%       (8)   MAPK model (from example_Benchmarking.m)

% Every combination of the models in 'testModels' and the reductions in
% 'reductionTypes' will be tested, and the results summarized in a table.
testModels = [1:8];
modelNames = {'2-Species Poisson','Repressilator','TV Repressilator',...
    'STL1','TXTL','Goutsias','Toggle2','MAPK'};

reductionTypes = {'POD2','Log Lump QSSA','QSSA'};    % {'Log Lump QSSA',
                           %  'Proper Orthogonal Decomposition','QSSA'};
reductionOrder = 30;
podTimeSetSize = 200;     % only needed for the POD reduction scheme
makePlots = true;        % plot full vs. reduced for every combination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose which type of model reduction to apply. Options include:
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
%   'Principle Orthogonal Decomposition' (POD) - 
%       The infintesimal generator is projected onto the orthonorml
%       basis spanned by a previous solution of the full master
%       equation at discrete time points.  This is done by perfoming
%       SVD on the solutions and choosing the output space
%       corresponding to the 'redOrder' largest singular values.
%       For best use, this reduction should be found using a
%       fine time resolution in the calculation of the FSP. The size of the
%       reduced model must be specified as 'reductionOrder'. Because the
%       POD requires a solution of the FSP, this reduction is usually only
%       helpful for situations where many solutions are needed (e.g.,
%       during model fitting).
%   'No Transform' - 
%       Test case where no reduction is applied.
%   Other Options that are only aplicable to time invariant systems:
%   'Log Lump QSSA' (LGQSSA) -- 
%       Same as LGSL, but where the distribution within each lump is
%       assumed to be distributed according to the quasi-steady state
%       assumption. 
%   'Dynamic Mode Decomposition' (DMD) -- 
%   (POD2) --
%       Extension of POD that also includes the time derivatives of the CME
%       (i.e., A*P(t_i)) solution in the set of vectors onto which the CME
%       is projected. Only the vectors corresponding top the top
%       'reductionOrder' singular values are kept for the projection.
%   'Eigen Decomposition'(ED) --
%       The infintesimal generator is projected onto the eigenvectors
%       corresponding to its least negative eigenvalues. The number of
%       eigenvectors for this projection is defined by 'reductionOrder'.      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Results are saved after each combination to a separate table for each
% model (named 'results_<modelName>'), using '-append' so that the results
% for one model can be updated without re-running the others.
resultsFile = 'example_SI_ModelReductionResults.mat';
iCombo = 0;

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
        case 2 % Constant 3-Species Repressilator Model
            Model1 = SSIT('Repressilator');
            Model1.tSpan = linspace(0,20,11);
            indTimes = [1,4,7,11];
            indSpecies = [3];
            qssaSpecies = [1,2];  % Species that are replaced with QSSA
        case 3 % Time-varying 3-Species Repressilator Model
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
        case 4 % Time-varying 4-state STL1 Model
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
            qssaSpecies = [5];  % Species that are replaced with QSSA

        case 5 % Time-varying Bursting TXTL Model
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
        case 6 % Goutsias Gene Regulation Model (from example_Benchmarking.m)
            Model1 = SSIT('Empty');
            Model1.parameters = {'k1',0.043;'k2',0.0007;'k3',0.0715;'k4',0.0039;'k5',0.012e9/(6.0221415e23*1e-15);...
                'k6',0.4791;'k7',0.00012e9/(6.0221415e23*1e-15);'k8',0.8765e-11;'k9',0.05e9/(6.0221415e23*1e-15);...
                'k10',0.5};
            Model1.species = {'D','DNA','M','RNA','DNAD','DNA2D'};
            Model1.stoichiometry = [0, 0,0, 0,-1, 1,-1, 1, 1, -1;...
                0, 0,0, 0,-1, 1, 0, 0, 0, 0;...
                1,-1,0, 0, 0, 0, 0, 0,-2, 2;...
                0, 0,1,-1, 0, 0, 0, 0, 0, 0;...
                0, 0,0, 0, 1,-1,-1, 1, 0, 0;...
                0, 0,0, 0, 0, 0, 1,-1, 0, 0];
            Model1.initialCondition = [6;2;2;0;0;0];
            % Reorder species DNA-RNA-Protein
            Model1.species = Model1.species([2,5,6,4,3,1]);
            Model1.initialCondition = Model1.initialCondition([2,5,6,4,3,1]);
            Model1.stoichiometry = Model1.stoichiometry([2,5,6,4,3,1],:);
            Model1.propensityFunctions = {'k1*RNA';'k2*M';'k3*DNAD';'k4*RNA';'k5*DNA*D';...
                'k6*DNAD';'k7*DNAD*D';'k8*DNA2D';'(k9/2)*M*(M-1)';'k10*D'};
            Model1 = Model1.addFSPConstraints(anticorrelatedPairs='all',correlatedPairs='all');
            Model1.tSpan = linspace(0,1000,11);
            indTimes = [1,4,7,11];
            indSpecies = [6];
            qssaSpecies = [1,2,3];  % Species that are replaced with QSSA
        case 7 % 2-Species Toggle Switch with Gene States (from example_Benchmarking.m)
            Model1 = SSIT('Empty');
            Model1.parameters = {'k1',40;'k2',20;'k3',1;'k4',1;'k5',1e-5;...
                'k6',3.5e-5;'k7',1;'k8',1};
            Model1.species = {'GeneA';'GeneB';'A';'B';'bGeneB';'bGeneA'};
            Model1.stoichiometry = [
                0,  0,  0,  0,  0, -1,  0,  1
                0,  0,  0,  0, -1,  0,  1,  0
                1,  0, -1,  0, -2,  0,  2,  0
                0,  1,  0, -1,  0, -2,  0,  2
                0,  0,  0,  0,  1,  0, -1,  0
                0,  0,  0,  0,  0,  1,  0, -1
                ];
            Model1.propensityFunctions = {'k1*GeneA';'k2*GeneB';'k3*A';'k4*B';...
                '(k5/2)*A*(A-1)*GeneB';'(k6/2)*B*(B-1)*GeneA';'k7*bGeneB';'k8*bGeneA'};
            Model1.initialCondition = [1;1;0;0;0;0];
            Model1 = Model1.addFSPConstraints(anticorrelatedPairs='all',correlatedPairs='all');
            Model1.tSpan = linspace(0,30,11);
            indTimes = [1,4,7,11];
            indSpecies = [3,4];
            qssaSpecies = [1,2,5,6];  % Species that are replaced with QSSA
        case 8 % MAPK Model (from example_Benchmarking.m)
            Model1 = SSIT();
            Model1.species = {'MKP3';'Kpp_MKP3';'KpT_MKP3_Y';'KpT_MKP3_T';...
                'K_MKP3_T';'KpY_MKP3';'K_MKP3_Y';'KpT_MEK';'KpT';'KpY_MEK';...
                'K_MEK_T';'KpY';'K_MEK_Y';'K';'MEK';'Kpp'};
            Model1.initialCondition = zeros(length(Model1.species),1);
            Model1.initialCondition(strcmp(Model1.species,'MKP3')) = 1;
            Model1.initialCondition(strcmp(Model1.species,'K')) = 3;
            Model1.parameters = {'s1',0.00024;'d1',0.0001;'s2',0.001;'d2',0.15;...
                's3',0.005;'k1',0.375;'k_1',1.0;'k2',0.06;'k3',0.375;'k_3',1.0;...
                'k4',4.5;'k5',0.375;'k_5',1.0;'k6',0.06;'k7',0.375;'k_7',1.0;...
                'k8',4.5;'h1',0.015;'h_1',1.0;'h2',0.032;'h3',0.31;'h_3',0.01;...
                'h4',0.01;'h_4',1.0;'h5',0.5;'h6',0.086;'h_6',0.0011;'h7',0.01;...
                'h_7',1.0;'h8',0.47;'h9',0.14;'h_9',0.0018};
            % Reactions: {propensity, {species, stoichiometry}}
            rxns = {'s2',{'MEK',1};
                'd2*MEK',{'MEK',-1};
                's3*Kpp',{'MEK',1};
                's1',{'K',1};
                'd1*K',{'K',-1};
                'd1*KpY',{'KpY',-1};
                'd1*KpT',{'KpT',-1};
                'd1*Kpp',{'Kpp',-1};
                'k1*K*MEK',{'K',-1;'MEK',-1;'K_MEK_Y',1};
                'k2*K_MEK_Y',{'K_MEK_Y',-1;'KpY',1;'MEK',1};
                'k3*KpY*MEK',{'KpY',-1;'MEK',-1;'KpY_MEK',1};
                'k_3*KpY_MEK',{'KpY_MEK',-1;'KpY',1;'MEK',1};
                'k4*KpY_MEK',{'KpY_MEK',-1;'Kpp',1;'MEK',1};
                'k5*K*MEK',{'K',-1;'MEK',-1;'K_MEK_T',1};
                'k_5*K_MEK_T',{'K_MEK_T',-1;'K',1;'MEK',1};
                'k6*K_MEK_T',{'K_MEK_T',-1;'KpT',1;'MEK',1};
                'k7*KpT*MEK',{'KpT',-1;'MEK',-1;'KpT_MEK',1};
                'k_7*KpT_MEK',{'KpT_MEK',-1;'KpT',1;'MEK',1};
                'k8*KpT_MEK',{'KpT_MEK',-1;'Kpp',1;'MEK',1};
                'h1*Kpp*MKP3',{'Kpp',-1;'MKP3',-1;'Kpp_MKP3',1};
                'h_1*Kpp_MKP3',{'Kpp_MKP3',-1;'Kpp',1;'MKP3',1};
                'h2*Kpp_MKP3',{'Kpp_MKP3',-1;'KpT_MKP3_Y',1};
                'h3*KpT_MKP3_Y',{'KpT_MKP3_Y',-1;'KpT',1;'MKP3',1};
                'h_3*KpT*MKP3',{'KpT',-1;'MKP3',-1;'KpT_MKP3_Y',1};
                'h4*KpT*MKP3',{'KpT',-1;'MKP3',-1;'KpT_MKP3_T',1};
                'h_4*KpT_MKP3_T',{'KpT_MKP3_T',-1;'KpT',1;'MKP3',1};
                'h5*KpT_MKP3_T',{'KpT_MKP3_T',-1;'K_MKP3_T',1};
                'h6*K_MKP3_T',{'K_MKP3_T',-1;'K',1;'MKP3',1};
                'h_6*K*MKP3',{'K',-1;'MKP3',-1;'K_MKP3_T',1};
                'h7*KpY*MKP3',{'KpY',-1;'MKP3',-1;'KpY_MKP3',1};
                'h_7*KpY_MKP3',{'KpY_MKP3',-1;'KpY',1;'MKP3',1};
                'h8*KpY_MKP3',{'KpY_MKP3',-1;'K_MKP3_Y',1};
                'h9*K_MKP3_Y',{'K_MKP3_Y',-1;'K',1;'MKP3',1};
                'h_9*K*MKP3',{'K',-1;'MKP3',-1;'K_MKP3_Y',1}};
            Model1.propensityFunctions = {};
            Model1.stoichiometry = [];
            for iRxn = 1:size(rxns,1)
                Model1 = Model1.addReaction(struct('propensity',rxns(iRxn,1),...
                    'stoichiometry',{rxns{iRxn,2}}));
            end
            % Add custom constraints based on sum of each monomer.
            Model1.customConstraintFuns = ...
                {'MKP3+Kpp_MKP3+KpT_MKP3_Y+KpT_MKP3_T+K_MKP3_T+KpY_MKP3+K_MKP3_Y';  % Sum of MKP3
                'K+K_MEK_Y+K_MEK_T+K_MKP3_T+K_MKP3_Y';  % Sum of K
                'Kpp+Kpp_MKP3'; % Sum of Kpp
                'KpY+KpY_MEK+KpY_MKP3'; % Sum of KpY
                'KpT+KpT_MEK+KpT_MKP3_Y+KpT_MKP3_T'; % Sum of KpT
                'K+K_MEK_Y+K_MEK_T+K_MKP3_T+K_MKP3_Y+KpY+KpY_MEK+KpY_MKP3+KpT+KpT_MEK+KpT_MKP3_Y+KpT_MKP3_T'}; % Sum of all
            Model1 = Model1.addFSPConstraints(anticorrelatedPairs='all',correlatedPairs='all');
            Model1.fspOptions.fspTol = 1e-3;
            Model1.tSpan = linspace(0,10,11);
            indTimes = [1,4,7,11];
            indSpecies = [16];
            qssaSpecies = [];  % Species that are replaced with QSSA
    end

    % Solve the original Model (once per model, for comparison)
    % Solve once to get the necessary FSP projection space.
    Model1 = Model1.formPropensitiesGeneral(['ModRedModel',num2str(testModel)]);
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
        iCombo = iCombo+1;
        status = 'OK';
        redModelSolveTime = NaN;
        maxKS = NaN;

        if strcmp(reductionType,'QSSA')&&isempty(qssaSpecies)
            status = 'QSSA Reduction Requires at least one QSSA Species';
        elseif strcmp(reductionType,'QSSA')&&length(Model1.species)==length(qssaSpecies)
            status = 'QSSA Reduction Requires at least one non-QSSA Species';
        else
            if ismember(reductionType,{'Proper Orthogonal Decomposition','POD','POD 2nd','POD2'})
                fspSoln2 = fspSolnFine;
            else
                fspSoln2 = fspSoln;
            end

            try
                % Solving the reduced models
                % Make a copy of the full model.
                Model2 = Model1;

                % Set the solver to use ModelReduction
                Model2.modelReductionOptions.useModReduction = true;
                % FSP expansion should be supressed when using Model Reductions

                % Set type and order of Model Recution
                Model2.modelReductionOptions.reductionType = reductionType;
                Model2.modelReductionOptions.reductionOrder = reductionOrder;
                Model2.modelReductionOptions.qssaSpecies = qssaSpecies;

                % Call SSIT to compute the model reduction transformation matrices.
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
            catch err
                status = err.message;
            end
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
            % Make Figures to compare the results. Here, we will plot the
            % original model in blue and the reduced model in red lines.
            figTitle = [modelNames{testModel},' - ',reductionType];

            % Means and standard deviations:
            f1 = figure(2*iCombo-1); clf;
            Model1.plotFSP(speciesNames=Model1.species(indSpecies),...
                plotType='meansAndDevs', lineProps={'linewidth',2},...
                Title=figTitle, TitleFontSize=26,...
                AxisLabelSize=20, TickLabelSize=20, YLabel='Molecule Count',...
                XLabel='Time', LegendFontSize=20, LegendLocation='southeast',...
                figureNums=f1);
            hold(gca,'on');
            Model2.plotFSP(speciesNames=Model2.species(indSpecies),...
                plotType='meansAndDevs', lineProps={'linewidth',2,'LineStyle','--'},...
                Title=figTitle, TitleFontSize=26,...
                AxisLabelSize=20, TickLabelSize=20, YLabel='Molecule Count',...
                XLabel='Time', LegendFontSize=20, LegendLocation='southeast',...
                figureNums=f1,Colors={'r'});

            % Marginal distributions:
            f2 = figure(2*iCombo); clf;
            Model1.plotFSP(speciesNames=Model1.species(indSpecies),...
                plotType='marginals', indTimes=indTimes,...
                lineProps={'linewidth',2}, Colors={'b'},figureNums=f2);
            hold(gca,'on');
            Model2.plotFSP(speciesNames=Model2.species(indSpecies),...
                plotType='marginals', indTimes=indTimes,...
                lineProps={'linewidth',2,'LineStyle','--'}, Colors={'r'},...
                figureNums=f2);
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
