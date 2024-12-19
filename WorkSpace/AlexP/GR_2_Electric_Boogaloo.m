%% Using the SSIT to fit Multiple Models and Data sets with Shared Parameters
% In this script, we show how multiple SSIT models and data sets can be fit
% simultaneously.  This is most useful in situations where:
%   1) the analysis considers different experimental conditions (e.g.,
%   different time points, different inducer concentrations, different
%   genetic mutations).
%   2) replica to replica variations are expected that would result in
%   slightly different parameter combinations

%% TEST NEW GR MODELS

close all 
clear
addpath(genpath('../../src'));
loadPrevious = false;
savedWorkspace = 'workspace_Boogaloo';
%addpath('tmpPropensityFunctions');

fitOptions = optimset('Display','iter','MaxIter',300);

GR = input('(0) GR-alpha only (base model);\n(0) PDO (GR-beta treated as distortion to GR-alpha data);\n(1) Convolve P(GR-alpha) and P(GR-beta);\n(1) The whole kit & caboodle (GR-alpha + GR-beta).\nChoose your destiny: ');

    % GR-alpha (base model) setup
        ModelGR_a = SSIT;
        ModelGR_a.species = {'cytGR_a';'nucGR_a'};
        ModelGR_a.initialCondition = [20;1];
        ModelGR_a.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'da1 * cytGR_a'; ...
                                         'ka1 * nucGR_a'; 'da2 * nucGR_a'};
        ModelGR_a.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        ModelGR_a.parameters = ({'MDex',5;'Dex0',100;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;...
                                 'ka1',0.01;'ba1',14e-5;'da1',1e-5;'da2',1e-6});
        ModelGR_a.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        ModelGR_a.summarizeModel

    % The log prior will be applied to the fit to multiple models as an additional constraint.
        log10PriorMean_a = [0.5 2 -3 -2 -2 -2 -4 -5 -6];
        log10PriorStd_a = 2*ones(1,9);

    % So it is left out of the prior, since we only want it to be calculated once.
        ModelGR_a.fittingOptions.logPrior = [];  
    
        ModelGR_a.fspOptions.initApproxSS = true;
        ModelGR_a.fittingOptions.modelVarsToFit = (3:9);
        
        ModelGR_a = ModelGR_a.formPropensitiesGeneral('EricRonModGR_a');
        ModelGR_a.customConstraintFuns = {'cytGR_a + nucGR_a'};
        % TODO - Alex - Constraints for removing stiff dimensions.
    
        [FSP_GR_a_Soln,ModelGR_a.fspOptions.bounds] = ModelGR_a.solve;
        [FSP_GR_a_Soln,ModelGR_a.fspOptions.bounds] = ModelGR_a.solve(FSP_GR_a_Soln.stateSpace);

    % STEP 0.B.2. -- Define GR parameters
        GR_a_pars = cell2mat(ModelGR_a.parameters(3:9,2))';  

    % STEP 0.B.3. -- Associate GR Data with Different Instances of Model (10,100nm Dex)
    GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
        '10','10',102,'GR Fit (10nM Dex)';...
        '100','100',103,'GR Fit (100nM Dex)'};
    ModelGRparameterMap = cell(1,size(GRfitCases,1));
    ModelGRfit = cell(1,size(GRfitCases,1));
    % ModelGRODEfit = cell(1,size(GRfitCases,1));
    for i=1:3
        ModelGRfit{i} = ModelGR_a.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
            {'nucGR_a','normgrnuc';'cytGR_a','normgrcyt'},...
            {'Dex_Conc',GRfitCases{i,2}});
        ModelGRfit{i}.parameters(9,:) = {'Dex0', str2num(GRfitCases{i,1})};
        ModelGRparameterMap(i) = {(1:7)};
        % parameters 3 - 9 refer to the parameter set that is relevant to
        % the entire class of models.  
    end
    
    % STEP 0.B.4. -- Make Guesses for the FSP bounds
    % This is sometimes necessary when using an uninduced steady state as the
    % initial condition. You need to guess a reasonable statespace or the
    % computation of the SS can be inaccurate.
    ModelGR.customConstraintFuns = {'cytGR+nucGR'};
    for i = 1:3
        boundGuesses{i} = [0;0;30;30;30];
        % First N are lower bounds.  Next N is upper bound.  Remaining are
        % custom.
    end

%% STEP 1 -- Fit GR Models.  
% STEP 1 will need to be rerun until satisfied.  Use fitMHiters as needed.
% TODO: Automate with statistics.
% Set for STEP1 -- Fit GR Models
fitIters = 30;
fitMHiters = 20;

for GR = 1:fitMHiters
    % STEP 1.A. -- Specify dataset time points.    
    for i = 1:3
        ModelGRfit{i}.tSpan = ModelGRfit{i}.dataSet.times;
    end

    % STEP 1.B. -- Specify log prior (NOTE: must transpose due to Matlab update that
    %     no longer correctly assumes format when adding single value vector to
    %     column vector).

    logPriorGR_a = @(x)-sum((log10(x)-log10PriorMean_a(1:7)').^2./(2*log10PriorStd_a(1:7)'.^2));

    % STEP 1.C. -- Combine all three GR models and fit using a single parameter set.
    for jj = 1:fitIters
        combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap,logPriorGR_a);
        combinedGRModel = combinedGRModel.initializeStateSpaces(boundGuesses);
        combinedGRModel = combinedGRModel.updateModels(GR_a_pars,false);
        GR_a_pars = combinedGRModel.maximizeLikelihood(GR_a_pars, fitOptions);
        save('EricModel_MMDex_GR2a','GR_a_pars') 
    end

    save('EricModelGR_MMDex_GR2a','GR_a_pars','combinedGRModel', 'ModelGRfit', 'log10PriorStd_a') 

    %% STEP 1.D. -- Compute FIM for GR parameters.
    combinedGRModel = combinedGRModel.computeFIMs([],'log');
    fimGR_withPrior = combinedGRModel.FIM.totalFIM+... % the FIM in log space.
        diag(1./(log10PriorStd_a(ModelGR_a.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.

    if min(eig(fimGR_withPrior))<1
        disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
        fimGR_withPrior = fimGR_withPrior + 1*eye(length(fimGR_withPrior));
    end
    covFree = fimGR_withPrior^-1;
    covFree = 0.5*(covFree+covFree');

    %% STEP 1.E. -- Run MH on GR Models.
    %GRpars = GRpars';
    MHFitOptions.thin=1;
    MHFitOptions.numberOfSamples=3000;
    MHFitOptions.burnIn=1000;
    MHFitOptions.progress=true;
    MHFitOptions.numChains = 1;

    % Use FIM computed above rather than making SSIT call 'useFIMforMetHast'
    % which forces SSIT.m to compute it within (no prior, etc.)
    MHFitOptions.useFIMforMetHast = false;
    MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);

    MHFitOptions.saveFile = 'TMPEricMHGR.mat';
    [~,~,MHResultsGR] = combinedGRModel.maximizeLikelihood(...
        GRpars, MHFitOptions, 'MetropolisHastings');
    %delete(MHFitOptions.saveFile)
    %MHResultsGR
    %%
    figNew = figure;
    ModelGR.plotMHResults(MHResultsGR,[],'log',[],figNew)
    for i = 1:7
        for j = i+1:7
            subplot(7,7,(i-1)*7+j-1)
            CH = get(gca,'Children');
            CH(1).Color=[1,0,1]; %
            CH(1).LineWidth = 3;
        end
    end
end 

switch GR
    case 0
        disp('You have chosen the base model, GR-alpha.')
    case 1
    % GR-beta setup
        disp('You have chosen to convolve P(GR-alpha) + P(GR-beta).')
        ModelGR_b = SSIT;
        ModelGR_b.species = {'cytGR_b';'nucGR_b'};
        ModelGR_b.initialCondition = [10;11];
        ModelGR_b.propensityFunctions = {'bb1'; 'db1 * cytGR_b'; 'kb1 * cytGR_b'; 'kb2 * nucGR_b'; 'db2 * nucGR_b'};
        ModelGR_b.stoichiometry = [1,-1,-1,1,0;...
                                   0,0,1,-1,-1];
        ModelGR_b.parameters = ({'kb1',0.01;'kb2',0.01;'bb1',14e-5;'db1',1e-5;'db2',1e-6});
        ModelGR_b.summarizeModel

        % The log prior will be applied to the fit to multiple models as an additional constraint.
        log10PriorMean_b = [-1 -1 -4 -5 -6];
        log10PriorStd_b = 2*ones(1,5);

        % So it is left out of the prior, since we only want it to be calculated once.
        ModelGR_b.fittingOptions.logPrior = [];  
    
        ModelGR_b.fspOptions.initApproxSS = true;
        ModelGR_b.fittingOptions.modelVarsToFit = (1:5);
        
        %ModelGR_b = ModelGR_b.formPropensitiesGeneral('EricRonModGR_b');
        ModelGR_b.customConstraintFuns = {'cytGR_b + nucGR_b'};
        % TODO - Alex - Constraints for removing stiff dimensions.
    
        [FSP_GR_b_Soln,ModelGR_b.fspOptions.bounds] = ModelGR_b.solve;
        [FSP_GR_b_Soln,ModelGR_b.fspOptions.bounds] = ModelGR_b.solve(FSP_GR_b_Soln.stateSpace);

    % STEP 0.B.2. -- Define GR parameters
        GR_b_pars = cell2mat(ModelGR_b.parameters(1:5,2))';
end
    
    

%%     STEP 1.F. -- Make Plots of GR Fit Results
makeGRPlots(combinedGRModel,GRpars)

save('EricModelGR_MMDex','GRpars','combinedGRModel','MHResultsGR') 
save('workspaceDec9_2024.mat','GRpars', 'ModelGRfit', 'combinedGRModel','MHResultsGR', 'log10PriorStd')