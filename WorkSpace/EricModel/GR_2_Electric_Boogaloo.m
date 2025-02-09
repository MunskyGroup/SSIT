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
addpath('tmpPropensityFunctions');
loadPrevious = false;
savedWorkspace = 'workspace_Boogaloo';

fitOptions = optimset('Display','iter','MaxIter',300);
 
%%
GR = input('(1) Convolve GR-alpha + GR-beta;\n(2) ODES for GR-alpha + GR-beta;\n(3) SSAs\nChoose your destiny: ');

switch GR
    case 1
    %% GR-alpha + GR-beta setup
        ModelGR_a = SSIT;
        ModelGR_b = SSIT;
        ModelGR_a.species = {'cytGR_a';'nucGR_a'};
        ModelGR_b.species = {'cytGR_b';'nucGR_b'};
        ModelGR_a.initialCondition = [5;1];
        ModelGR_b.initialCondition = [5;1];
        ModelGR_a.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'dc * cytGR_a'; 'kn2c * nucGR_a'; 'dn * nucGR_a'};
        ModelGR_b.propensityFunctions = {'kb1 * cytGR_b'; 'bb1'; 'dc * cytGR_b'; 'kn2c * nucGR_b'; 'dn * nucGR_b'};
        ModelGR_a.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        ModelGR_b.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        ModelGR_a.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'ba1',14e-5;...
                                'MDex',5;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;'Dex0',100});
        ModelGR_b.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'kb1',0.01;'bb1',14e-5});
        ModelGR_a.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        disp('The GR-alpha model is: ')
        ModelGR_a.summarizeModel
        disp('and the GR-beta model is: ')
        ModelGR_b.summarizeModel

        ModelGR_a.fspOptions.initApproxSS = true;
        ModelGR_b.fspOptions.initApproxSS = true;

        ModelGR_a.fittingOptions.modelVarsToFit = (1:8);
        ModelGR_b.fittingOptions.modelVarsToFit = (1:5);
         
        ModelGR_a = ModelGR_a.formPropensitiesGeneral('GR_a');
        ModelGR_b = ModelGR_b.formPropensitiesGeneral('GR_b');

        ModelGR_a.customConstraintFuns = {'cytGR_a+nucGR_a'};
        ModelGR_b.customConstraintFuns = {'cytGR_b+nucGR_b'};
        % % TODO - Alex - Constraints for removing stiff dimensions.

        %% Compute FSP solutions
        [FSPsoln_a,ModelGR_a.fspOptions.bounds] = ModelGR_a.solve;
        [FSPsoln_a,ModelGR_a.fspOptions.bounds] = ModelGR_a.solve(FSPsoln_a.stateSpace);

        [FSPsoln_b,ModelGR_b.fspOptions.bounds] = ModelGR_b.solve;
        [FSPsoln_b,ModelGR_b.fspOptions.bounds] = ModelGR_b.solve(FSPsoln_b.stateSpace);
        
        %% Convolution
        for f=1:(max(numel(FSPsoln_a.fsp),numel(FSPsoln_b.fsp)))
            % f is time point, so solution tensors are FSP probabilities across states for each time point
            conv2solnTensor{f} = conv2(double(FSPsoln_a.fsp{f}.p.data),double(FSPsoln_b.fsp{f}.p.data));
            figure(f)
            contourf(log10(conv2solnTensor{f}))
            hold on
        end
        %% check
        for g=1:f % f is number of time points 
            fspsoln_sptensor_a{g} = double(FSPsoln_a.fsp{g}.p.data);
            fspsoln_sptensor_b{g} = double(FSPsoln_b.fsp{g}.p.data);
            figure(g)
            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a{g}))
            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b{g}))
            subplot(1,3,3)
            contourf(log10(conv2solnTensor{g}))
            hold on
        end
         
        %% TODO: clean up scattered data on contour plots

        data = readtable('../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv');

        % Filter the data if desired for particular Time_index and/or Dex_Conc
        filteredData = data(data.time == 180 & data.Dex_Conc == 100, :);

        % Randomly sample 1000 points from the filtered data
        numSamples = min(1000, height(filteredData)); % Ensure we don’t sample more than available data
        randomIndices = randperm(height(filteredData), numSamples);
        normgrnuc = filteredData.normgrnuc(randomIndices);
        normgrcyt = filteredData.normgrcyt(randomIndices);

        % Plot contourf figures and add scatter points
        for g=1:f
            fspsoln_sptensor_a{g} = double(FSPsoln_a.fsp{g}.p.data);
            fspsoln_sptensor_b{g} = double(FSPsoln_b.fsp{g}.p.data);
    
            figure(g)
    
            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a{g}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b{g}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,3)
            contourf(log10(conv2solnTensor{g}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
        end
                

        %% Fit GR_a to 
        %% Define GR parameters
        GRpars_a = cell2mat(ModelGR_a.parameters(1:8,2))'; 
        GRpars_b = cell2mat(ModelGR_b.parameters(1:5,2))';


        % The log prior will be applied to the fit to multiple models as an additional constraint.
        log10PriorMean_a = [-2 -5 -6 -5 0.5 -3 -3 -2];
        log10PriorMean_b = [-2 -5 -6 -2 -5];
        log10PriorStd_a = 2*ones(1,8);
        log10PriorStd_b = 2*ones(1,5);
         
        % So it is left out of the prior, since we only want it to be calculated once.
        ModelGR_a.fittingOptions.logPrior = [];  
        ModelGR_b.fittingOptions.logPrior = [];          
    
        %% Associate GR Data with Different Instances of Model (1,10,100nm Dex)
        GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
            '10','10',102,'GR Fit (10nM Dex)';...
            '100','100',103,'GR Fit (100nM Dex)'};

        ModelGRparameterMap_a = cell(1,size(GRfitCases,1));
        ModelGRfit_a = cell(1,size(GRfitCases,1));
        ModelGRfit_b = cell(1,1);

        % Load data for GR-beta
        ModelGRfit_b{1} = ModelGR_b.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR_b','normgrnuc';'cytGR_b','normgrcyt'},{'Dex_Conc','100'});

        % ModelGRODEfit = cell(1,size(GRfitCases,1));
        
        % Load data, fit 3 Dex concentrations to GR-alpha's 'Dex0' parameter,
        % and set solution scheme to FSP.
        for i=1:3
            ModelGRfit_a{i} = ModelGR_a.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR_a','normgrnuc';'cytGR_a','normgrcyt'},...
                {'Dex_Conc',GRfitCases{i,2}});
            ModelGRfit_a{i}.parameters(9,:) = {'Dex0', str2num(GRfitCases{i,1})};
            ModelGRparameterMap_a(i) = {(1:8)}; 
        end
 
        % Make Guesses for the FSP bounds
        % This is sometimes necessary when using an uninduced steady state as the
        % initial condition. You need to guess a reasonable statespace or the
        % computation of the SS can be inaccurate.
        for i = 1:3
            boundGuesses{i} = [0;0;30;30;30];
            % First N are lower bounds.  Next N is upper bound.  Remaining are
            % custom.
        end

    % Specify log prior (NOTE: must transpose due to Matlab update that
    %     no longer correctly assumes format when adding single value vector to
    %     column vector).

    logPriorGR_a = @(x)-sum((log10(x)-log10PriorMean_a(1:8)').^2./(2*log10PriorStd_a(1:8)'.^2));
    logPriorGR_b = @(x)-sum((log10(x)-log10PriorMean_b(1:5)').^2./(2*log10PriorStd_b(1:5)'.^2));
    
    
    %% Fit GR Models.  
    % This step will need to be rerun until satisfied.  Use fitMHiters as needed.
    % TODO: Automate with statistics.

    fitIters = 3; % 3 concentrations of Dex (1,10,100nM)
    fitMHiters = 2; % Adjust as needed
    
    for GR = 1:fitMHiters
        % Specify dataset time points for GR-alpha.    
        for i = 1:3
            ModelGRfit_a{i}.tSpan = ModelGRfit_a{i}.dataSet.times;
        end

        % Specify dataset time points for GR-beta.
        ModelGRfit_b{1}.tSpan = ModelGRfit_b{1}.dataSet.times;

        %% Solve for combined GR-alpha model for three Dex_Conc=Dex0 and fit using a single parameter set.
        for jj = 1:fitIters
            % Solve
            combinedGRModel_a = SSITMultiModel(ModelGRfit_a,ModelGRparameterMap_a,logPriorGR_a);
            combinedGRModel_a = combinedGRModel_a.initializeStateSpaces(boundGuesses);
            combinedGRModel_a = combinedGRModel_a.updateModels(GRpars_a,false);
            GRpars_a = combinedGRModel_a.maximizeLikelihood(GRpars_a, fitOptions);
            save('combinedGRModel_a','GRpars_a') 
        end
        
        %% Solve for GR-beta model.
        ModelGRfit_b{1}.fspOptions.fspTol = 1e-4;
        ModelGRfit_b{1}.fspOptions.bounds = boundGuesses{1};
        [GR_b_fspSoln,ModelGRfit_b{1}.fspOptions.bounds] = ModelGRfit_b{1}.solve;
        [GR_b_fspSoln,ModelGRfit_b{1}.fspOptions.bounds] = ModelGRfit_b{1}.solve(GR_b_fspSoln.stateSpace);
        GRpars_b = ModelGRfit_b{1}.maximizeLikelihood(GRpars_b, fitOptions);
        save('ModelGRfit_b','GRpars_b','GR_b_fspSoln')             

        save('GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln')
        

        %% Compute FIM 
        % for ModelGR_a parameters.
        combinedGRModel_a = combinedGRModel_a.computeFIMs([],'log');

        fimGR_a_withPrior = combinedGRModel_a.FIM.totalFIM+... % the FIM in log space.
            diag(1./(log10PriorStd_a(ModelGR_a.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.

        % for ModelGR_b parameters.
        %ModelGR_b_fimResults = ModelGRfit_b{1}.computeFIM([],'log'); % Compute individual FIMs

        %fimGR_b_withPrior = ModelGRfit_b{1}.totalFim(ModelGR_b_fimResults,ModelGRfit_b{1}.dataSet.nCells,diag(1./(log10PriorStd_b(ModelGRfit_b{1}.fittingOptions.modelVarsToFit)*log(10)).^2));  % Add prior in log space.

        %
        if min(eig(fimGR_a_withPrior))<1
            disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
            fimGR_a_withPrior = fimGR_a_withPrior + 1*eye(length(fimGR_a_withPrior));
        end
        fimGR_a_covFree = fimGR_a_withPrior^-1;
        fimGR_a_covFree = 0.5*(fimGR_a_covFree+fimGR_a_covFree');

        % if min(eig(fimGR_b_withPrior{1}))<1
        %     disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
        %     fimGR_b_withPrior{1} = fimGR_b_withPrior{1} + 1*eye(length(fimGR_b_withPrior{1}));
        % end
        % fimGR_b_covFree = fimGR_b_withPrior{1}^-1;
        % fimGR_b_covFree = 0.5*(fimGR_b_covFree+fimGR_b_covFree');

    
        %% Run MH on GR Models.
        %GRpars = GRpars';
        % GR_alpha
        MHfitOptions_a.thin=1;
        MHfitOptions_a.numberOfSamples=1000;
        MHfitOptions_a.burnIn=100;
        MHfitOptions_a.progress=true;
        MHfitOptions_a.numChains = 1;
    
        % Use FIM computed above rather than making SSIT call 'useFIMforMetHast'
        % which forces SSIT.m to compute it within (no prior, etc.)
        MHfitOptions_a.useFIMforMetHast = false;
        MHfitOptions_a.proposalDistribution=@(x)mvnrnd(x,fimGR_a_covFree);
    
        MHfitOptions_a.saveFile = 'TMPEricMHGR_a.mat';
        [~,~,MHResultsGR_a] = combinedGRModel_a.maximizeLikelihood(...
            GRpars_a, MHfitOptions_a, 'MetropolisHastings');
        %delete(MHFitOptions.saveFile)
        %MHResultsGR
        %
        figNew = figure;
        ModelGR_a.plotMHResults(MHResultsGR_a,[],'log',[],figNew)
        for i = 1:7
            for j = i+1:7
                subplot(7,7,(i-1)*7+j-1)
                CH = get(gca,'Children');
                CH(1).Color=[1,0,1]; %
                CH(1).LineWidth = 3;
            end
        end
        %% GR_beta
        % MHFitOptions_b.thin=1;
        % MHFitOptions_b.numberOfSamples=1000;
        % MHFitOptions_b.burnIn=100;
        % MHFitOptions_b.progress=true;
        % MHFitOptions_b.numChains = 1;
        % 
        % % Use FIM computed above rather than making SSIT call 'useFIMforMetHast'
        % % which forces SSIT.m to compute it within (no prior, etc.)
        % MHFitOptions_b.useFIMforMetHast = false;
        % MHFitOptions_b.proposalDistribution=@(y)mvnrnd(y,fimGR_b_covFree);
        % 
        % MHFitOptions_b.saveFile = 'TMPEricMHGR_b.mat';
        % [~,~,MHResultsGR_b] = ModelGRfit_b{1}.maximizeLikelihood(...
        %     GRpars_b, MHFitOptions_b, 'MetropolisHastings');
        % %delete(MHFitOptions.saveFile)
        % %MHResultsGR
        % %
        % figNew = figure;
        % ModelGR_b.plotMHResults(MHResultsGR_b,[],'log',[],figNew)
        % for i = 1:7
        %     for j = i+1:7
        %         subplot(7,7,(i-1)*7+j-1)
        %         CH = get(gca,'Children');
        %         CH(1).Color=[1,0,1]; %
        %         CH(1).LineWidth = 3;
        %     end
        % end
    end

    %% Convolution, contour plots, and scattered data
    load('combinedGR_a_fspSolns.mat')
    data = readtable('../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv'); 

    tspan = [0,10,30,50,75,150,180];

    for model=1:size(GRfitCases,1)
        dataDex = data(data.Dex_Conc == str2double(GRfitCases{model}), :);
        for t=1:max(numel(fspSolnsSMM(model).fsp),numel(GR_b_fspSoln.fsp))
            % Figure index
            figIndex = (model - 1) * 7 + t;
            % Filter data           
            dataDexTime = dataDex(dataDex.time == tspan(t),:);

            % Randomly sample 1000 points from the filtered data
            numSamples = min(1000, height(dataDexTime)); % Ensure we don’t sample more than available data
            % check
            if numSamples==0
                GRfitCases{model}
                tspan(t)
            end
            randomIndices = randperm(height(dataDexTime), numSamples);
            normgrnuc = dataDexTime.normgrnuc(randomIndices);
            normgrcyt = dataDexTime.normgrcyt(randomIndices);
                
            % t is time point, so solution tensors are FSP probabilities across states for each time point
            conv2solnTensor_postData{t} = conv2(double(fspSolnsSMM(model).fsp{t}.p.data),double(GR_b_fspSoln.fsp{t}.p.data));
            figure(figIndex)
            contourf(log10(conv2solnTensor_postData{t}))
            hold on 
            scatter(normgrcyt, normgrnuc, 'r', 'filled')

            % Add subplots (can be commented out)
            fspsoln_sptensor_a_postData{t} = double(fspSolnsSMM(model).fsp{t}.p.data);
            fspsoln_sptensor_b_postData{t} = double(GR_b_fspSoln.fsp{t}.p.data);

            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a_postData{t}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b_postData{t}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,3)
            contourf(log10(conv2solnTensor_postData{t}))
            hold on
            scatter(normgrcyt, normgrnuc, 'r', 'filled') % Red dots
        end
    end

    save('conv2solnTensor_postData','GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln', 'fspSolnsSMM')

    dusp1 = input('(0) GR-beta turns off the DUSP1 gene;\n(1) GR-beta has no effect on DUSP1;\nChoose your destiny: ');

    switch dusp1
        case 0 
            %% Extend model for DUSP1
            % GR-beta turns off the DUSP1 gene
            ModelGRDusp100nM = ModelGRfit_a{3};
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'offGene'},2);
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'onGene'},0);
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'rna'},5);
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'cytGR_b'},5);
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'nucGR_b'},1);

            ModelGRDusp100nM.propensityFunctions{6,1} = 'kon*offGene*nucGR_a';
            ModelGRDusp100nM.propensityFunctions{7,1} = 'koff*onGene*nucGR_b';
            ModelGRDusp100nM.propensityFunctions{8,1} = 'kr*onGene';
            ModelGRDusp100nM.propensityFunctions{9,1} = 'dr*rna';
            ModelGRDusp100nM.propensityFunctions{10,1} = 'kb1 * cytGR_b';
            ModelGRDusp100nM.propensityFunctions{11,1} = 'bb1'; 
            ModelGRDusp100nM.propensityFunctions{12,1} = 'dc * cytGR_b';
            ModelGRDusp100nM.propensityFunctions{13,1} = 'kn2c * nucGR_b';
            ModelGRDusp100nM.propensityFunctions{14,1} = 'dn * nucGR_b';

            ModelGRDusp100nM.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0,0,0,0,0;...
                                               1,0,0,-1,-1,0,0,0,0,0,0,0,0,0;...
                                               0,0,0,0,0,-1,1,0,0,0,0,0,0,0;...
                                               0,0,0,0,0,1,-1,0,0,0,0,0,0,0;...
                                               0,0,0,0,0,0,0,1,-1,0,0,0,0,0;...
                                               0,0,0,0,0,0,0,0,0,-1,1,-1,1,0;...
                                               0,0,0,0,0,0,0,0,0,1,0,0,-1,-1];

            ModelGRDusp100nM.parameters(10,:) = {'koff',0.1};
            ModelGRDusp100nM.parameters(11,:) = {'kon',0.1};
            ModelGRDusp100nM.parameters(12,:) = {'kr',1};
            ModelGRDusp100nM.parameters(13,:) = {'dr',0.02};

            ModelGRDusp100nM.parameters(14,:) = {'kb1',GRpars_b(4)};
            ModelGRDusp100nM.parameters(15,:) = {'bb1',GRpars_b(5)};

            ModelGRDusp100nM.useHybrid = true;
            ModelGRDusp100nM.hybridOptions.upstreamODEs = {'cytGR_a','nucGR_a','cytGR_b','nucGR_b'};

            ModelGRDusp100nM.solutionScheme = 'FSP';
            ModelGRDusp100nM.customConstraintFuns = [];
            %ModelGRDusp100nM.fspOptions.bounds = [0;0;0;2;2;400];
            ModelGRDusp100nM.fittingOptions.modelVarsToFit = 10:13;
            ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('EricModDusp1_1');
            log10PriorMean = [-2 -5 -6 -5 0.5 -3 -3 -2 2 -1 -1 0 -2 -2 -5];
            log10PriorStd = 2*ones(1,15);
            duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(9:12)).^2./(2*log10PriorStd(9:12).^2));
            ModelGRDusp100nM.fittingOptions.logPrior = duspLogPrior;

            ModelGRDusp100nM.summarizeModel

        %% 
            Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
            ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
            ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
            for i = 1:size(Dusp1FitCases,1)
                 ModelDusp1Fit{i} = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
                     {'rna','totalNucRNA'},...
                     {'Dex_Conc','100'});
                 ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
                 ModelDusp1parameterMap{i} = (1:4);
                 % Set Dex concentration.
                 ModelDusp1Fit{i}.parameters{9,2} = str2num(Dusp1FitCases{i,1});
                 ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
            end
            DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];

        %% Load data
            ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
                                                            {'rna','RNA_DUSP1_nuc'},{'Dex_Conc','100'});
            DUSP1pars = [ModelGRDusp100nM.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
    
        %% Fit DUSP1 model at 100nM Dex.  
        % This will likely need to be rerun after Metropolis-Hastings Search .
        % Use fitMHiters to run until satisfied. TODO: Automate by computing statistics.
        fitMHiters = 2;
        
        for DS = 1:fitMHiters
            for i = 1:fitIters
                fitOptions.suppressFSPExpansion = true;
                DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
                    DUSP1pars, fitOptions);
                ModelGRDusp100nM.parameters(10:13,2) = num2cell(DUSP1pars);
                save('EricModelDusp1_MMDex','GRpars_a','GRpars_b','DUSP1pars') 
            end
        
            %%  Plot predictions for other Dex concentrations.
            %showCases = [1,1,1,1];
            %makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1pars,Dusp1FitCases,showCases)
        
            %% Sample uncertainty for Dusp1 Parameters
            %   Compute sensitivity of the FSP solution
            ModelGRDusp100nM.solutionScheme = 'fspSens';
            sensSoln = ModelGRDusp100nM.solve();
            ModelGRDusp100nM.solutionScheme = 'FSP';
            %   Compute FIM
            %       define which species in model are not observed.
            ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};
            % compute the FIM
            fimResults = ModelGRDusp100nM.computeFIM(sensSoln.sens,'log');
            % In the following, the log-prior is used as a prior co-variance matrix.
            % This will be used in the FIM calculation as an FIM without new evidence 
            % being set equal to the inverse of this covariance matrix.  More rigorous
            %% justification is needed to support this heuristic.
            fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
                diag(log10PriorStd(1:15).^2));
            FIMfree = fimTotal{1}(1:4,1:4);
            if min(eig(FIMfree))<1
                disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
                FIMfree = FIMfree + 1*eye(length(FIMfree));
            end
            covFree = FIMfree^-1;
            covFree = 0.5*(covFree+covFree');
        
            %% Run Metropolis Hastings Search
            if loadPrevious
                MHDusp1File = 'MHDusp1_Dec92024';
                load(MHDusp1File)
            else
                MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
                MHFitOptions.thin=1;
                MHFitOptions.numberOfSamples=2000;
                MHFitOptions.burnIn=0;
                MHFitOptions.progress=true;
                MHFitOptions.numChains = 1;
                MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
                [DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
                    [], MHFitOptions, 'MetropolisHastings');
                delete('TMPEricMHDusp1.mat')
                ModelGRDusp100nM.parameters(1:4,2) = num2cell(DUSP1pars);
            end
            
            save('workspace_Feb4_2024.mat','ModelGRDusp100nM','DUSP1pars','fimTotal','sensSoln','combinedGRModel_a','MHResultsDusp1')
        
            %% Plot the MH results
            figNew = figure;
            ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
            for j = 1:3
                for k = i:3
                    subplot(3,3,(j-1)*3+k)
                    CH = get(gca,'Children');
                    CH(1).Color=[1,0,1]; %
                    CH(1).LineWidth = 3;
                end
            end
        end

        case 1  %% GR-beta has no effect on DUSP1

        %% Add DUSP1 to the existing GR model.
            % Copy parameters from the 100nM Dex stim case in GR.
            fitIters = 3;
            ModelGRDusp100nM = ModelGRfit_a{3};
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'offGene'},2);
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'onGene'},0);
            ModelGRDusp100nM = ModelGRDusp100nM.addSpecies({'rna'},5);
            ModelGRDusp100nM.propensityFunctions{6,1} = 'kon*offGene*nucGR_a';
            ModelGRDusp100nM.propensityFunctions{7,1} = 'koff*onGene';
            ModelGRDusp100nM.propensityFunctions{8,1} = 'kr*onGene';
            ModelGRDusp100nM.propensityFunctions{9,1} = 'dr*rna';
            ModelGRDusp100nM.stoichiometry = [-1,1,-1,1,0,0,0,0,0;...
                                        1,0,0,-1,-1,0,0,0,0;...
                                        0,0,0,0,0,-1,1,0,0;...
                                        0,0,0,0,0,1,-1,0,0;...
                                        0,0,0,0,0,0,0,1,-1];
            ModelGRDusp100nM.parameters(10,:) = {'koff',0.1};
            ModelGRDusp100nM.parameters(11,:) = {'kon',0.1};
            ModelGRDusp100nM.parameters(12,:) = {'kr',1};
            ModelGRDusp100nM.parameters(13,:) = {'dr',0.02};
            ModelGRDusp100nM.useHybrid = true;
            ModelGRDusp100nM.hybridOptions.upstreamODEs = {'cytGR_a','nucGR_a'};
            ModelGRDusp100nM.solutionScheme = 'FSP';
            ModelGRDusp100nM.customConstraintFuns = [];
            ModelGRDusp100nM.fspOptions.bounds = [0;0;0;2;2;400];
            ModelGRDusp100nM.fittingOptions.modelVarsToFit = 10:13;
            ModelGRDusp100nM = ModelGRDusp100nM.formPropensitiesGeneral('EricModDusp1');
            log10PriorMean = [-2 -5 -6 -5 0.5 -3 -3 -2 2 -1 -1 0 -2];
            log10PriorStd = 2*ones(1,13);
            duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(10:13)).^2./(2*log10PriorStd(10:13).^2));
            ModelGRDusp100nM.fittingOptions.logPrior = duspLogPrior;
        
        %% 
            Dusp1FitCases = {'100','100',201,'DUSP1 Fit (100nM Dex)'};
            ModelDusp1Fit = cell(size(Dusp1FitCases,1),1);
            ModelDusp1parameterMap = cell(1,size(GRfitCases,1));
            for i = 1:size(Dusp1FitCases,1)
                 ModelDusp1Fit{i} = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
                     {'rna','totalNucRNA'},...
                     {'Dex_Conc','100'});
                 ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
                 ModelDusp1parameterMap{i} = (1:4);
                 % Set Dex concentration.
                 ModelDusp1Fit{i}.parameters{13,2} = str2num(Dusp1FitCases{i,1});
                 ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
            end
            DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];

        %%
        ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
            {'rna','totalNucRNA'},{'Dex_Conc','100'});
        DUSP1pars = [ModelGRDusp100nM.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];
    
        %% Fit DUSP1 model at 100nM Dex.  
        % This will likely need to be rerun after Metropolis-Hastings Search .
        % Use fitMHiters to run until satisfied. TODO: Automate by computing statistics.
        fitMHiters = 2;
        
        for DS = 1:fitMHiters
            for i = 1:fitIters
                fitOptions.suppressFSPExpansion = true;
                DUSP1pars = ModelGRDusp100nM.maximizeLikelihood(...
                    DUSP1pars, fitOptions);
                ModelGRDusp100nM.parameters(10:13,2) = num2cell(DUSP1pars);
                save('EricModelDusp1_MMDex','GRpars_a','GRpars_b','DUSP1pars') 
            end
        
            %% Plot predictions for other Dex concentrations.
            %showCases = [1,1,1,1];
            %makePlotsDUSP1({ModelGRDusp100nM},ModelGRDusp100nM,DUSP1pars,Dusp1FitCases,showCases)
        
            %% Sample uncertainty for Dusp1 Parameters
            %  Compute sensitivity of the FSP solution
            ModelGRDusp100nM.solutionScheme = 'fspSens';
            sensSoln = ModelGRDusp100nM.solve();
            ModelGRDusp100nM.solutionScheme = 'FSP';
            %  Compute FIM
            %       define which species in model are not observed.
            ModelGRDusp100nM.pdoOptions.unobservedSpecies = {'offGene';'onGene'};
            % compute the FIM
            fimResults = ModelGRDusp100nM.computeFIM(sensSoln.sens,'log');
            % In the following, the log-prior is used as a prior co-variance matrix.
            % This will be used in the FIM calculation as an FIM without new evidence 
            % being set equal to the inverse of this covariance matrix.  More rigorous
            %% justification is needed to support this heuristic.
            fimTotal = ModelGRDusp100nM.evaluateExperiment(fimResults,ModelGRDusp100nM.dataSet.nCells,...
                diag(log10PriorStd(1:13).^2));
            FIMfree = fimTotal{1}(10:13,10:13);
            if min(eig(FIMfree))<1
                disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
                FIMfree = FIMfree + 1*eye(length(FIMfree));
            end
            covFree = FIMfree^-1;
            covFree = 0.5*(covFree+covFree');
        
            %%      STEP 2.D.3. -- Run Metropolis Hastings Search
            if loadPrevious
                MHDusp1File = 'MHDusp1_Feb_2025';
                load(MHDusp1File)
            else
                MHFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
                MHFitOptions.thin=1;
                MHFitOptions.numberOfSamples=2000;
                MHFitOptions.burnIn=0;
                MHFitOptions.progress=true;
                MHFitOptions.numChains = 1;
                MHFitOptions.saveFile = 'TMPEricMHDusp1.mat';
                [DUSP1pars,~,MHResultsDusp1] = ModelGRDusp100nM.maximizeLikelihood(...
                    [], MHFitOptions, 'MetropolisHastings');
                delete('TMPEricMHDusp1.mat')
                ModelGRDusp100nM.parameters(10:13,2) = num2cell(DUSP1pars);
            end
            
            save('workspace_Feb4_2024.mat','ModelGRDusp100nM','DUSP1pars','fimTotal','sensSoln','combinedGRModel_a','ModelGRfit_a','ModelGRfit_b','MHResultsDusp1', 'conv2solnTensor_postData')
        
            %% Plot the MH results
            figNew = figure;
            ModelGRDusp100nM.plotMHResults(MHResultsDusp1,[],'log',[],figNew)
            for j = 1:3
                for k = i:3
                    subplot(3,3,(j-1)*3+k)
                    CH = get(gca,'Children');
                    CH(1).Color=[1,0,1]; %
                    CH(1).LineWidth = 3;
                end
            end
        end
        
        %% Save results
        varNames = unique({'ModelGR_a'
            'ModelGR_b'
            'GRfitCases'
            'log10PriorMean'
            'log10PriorStd'
            'GRpars_a'
            'GRpars_b'
            'ModelGRparameterMap'
            'combinedGR_a'
            'ModelGRDusp100nM'
            'duspLogPrior'
            'DUSP1pars'
            'ModelGRfit_a'
            'ModelGRfit_b'
            'fimResults'
            'MHResultsDusp1'
            'ModelGRDusp100nM'
            });
        
        save('workspace_Feb_2025',varNames{:}) % WARNING: THIS OVERWRITE THE PREVIOUSLY SAVED WORKSPACE - TODO: FIX
    end

    %save('conv2solnTensor_postData','GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln', 'fspSolnsSMM',...
    %    'fimGR_a_withPrior','fimGR_b_withPrior','ModelGR_b_fimResults','fimGR_a_covFree','fimGR_b_covFree','MHResultsGR_a','MHFitOptions','ModelGRfit');
    save('conv2solnTensor_postData','GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln', 'fspSolnsSMM',...
        'fimGR_a_withPrior','fimGR_a_covFree','MHResultsGR_a','MHFitOptions','ModelGRfit');

    %% Create new objective function combining all of the previous ones.
     % Remove all priors from individual models.
    for i=1:3
        ModelGRfit_a{i}.fittingOptions.logPrior = [];
        ModelGRfit_a{i}.tSpan = ModelGRfit_a{i}.dataSet.times;
    end
    ModelGRfit_b{1}.fittingOptions.logPrior = [];
    ModelGRfit_b{1}.tSpan = ModelGRfit_b{1}.dataSet.times;
    ModelGRDusp100nM.fittingOptions.logPrior = [];
    fullPars = [ModelGRDusp100nM.parameters{1:15,2}];

    %% Fit all objective functions at once.
    log10PriorMean = [-2 -5 -6 ...          % GR pars
                      -5 0.5 -3 -3 -2 ...   % GR-alpha pars
                       2 ...              % Dex
                      -1 -1 0 -2 ...        % Dusp1 pars
                      -2 -5];               % GR-beta pars
    log10PriorStd = 2*ones(1,15);
    %duspLogPrior = @(x)-sum((log10(x(:))'-log10PriorMean(9:12)).^2./(2*log10PriorStd(9:12).^2));
    logPriorAll = @(x)-sum((log10(x)-log10PriorMean(1:15)).^2./(2*log10PriorStd(1:15).^2));

        % extendedMod.fittingOptions.modelVarsToFit = [1:12,14,15];
        Organization = {ModelGRfit_a{3},[1:8],[1:8],'computeLikelihood',1;...  %% TODO: get convol. FSP solns to feed in
            ModelGRfit_b{1},[1:5],[1:5],'computeLikelihood',1;...
            ModelGRDusp100nM,[1:8,10:15],[1:8,10:15],'computeLikelihood',1};
            %extendedMod,[1:12,14:15],[1:14],'computeLikelihoodODE',0.01};
        Organization = getTotalFitErr(Organization,fullPars,true);
        getTotalFitErr(Organization,fullPars,false)
        objAll = @(x)-getTotalFitErr(Organization,exp(x),false)-logPriorAll(exp(x));
        for jj = 1:fitIters
            fullParsLog = log(fullPars);
            fullPars = exp(fminsearch(objAll,fullParsLog,fitOptions));
        end    
    
    %% compute likelihood
    %    logLikelihood = sum(log(conv2solnTensor{f}(PAs>0)).*PAs((PAs>0)),"all");
    %logLikelihood = sum(log(conv2solnTensor_postData{t}(fspsoln_sptensor_a_postData{t}>0)).*fspsoln_sptensor_a_postData{t}((fspsoln_sptensor_a_postData{t}>0)),"all");

    case 2 
        %% ODEs        
        ModelGR_a = SSIT;
        ModelGR_b = SSIT;
        ModelGR_a.species = {'cytGR_a';'nucGR_a'};
        ModelGR_b.species = {'cytGR_b';'nucGR_b'};
        ModelGR_a.initialCondition = [5;1];
        ModelGR_b.initialCondition = [5;1];
        ModelGR_a.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'dc * cytGR_a'; 'kn2c * nucGR_a'; 'dn * nucGR_a'};
        ModelGR_b.propensityFunctions = {'kb1 * cytGR_b'; 'bb1'; 'dc * cytGR_b'; 'kn2c * nucGR_b'; 'dn * nucGR_b'};
        ModelGR_a.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        ModelGR_b.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        ModelGR_a.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'ba1',14e-5;...
                                'MDex',5;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;'Dex0',100});
        ModelGR_b.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'kb1',0.01;'bb1',14e-5});
        ModelGR_a.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        disp('The GR-alpha model is: ')
        ModelGR_a.summarizeModel
        disp('and the GR-beta model is: ')
        ModelGR_b.summarizeModel

        ModelGR_a.fspOptions.initApproxSS = true;
        ModelGR_b.fspOptions.initApproxSS = true;

        ModelGR_a.fittingOptions.modelVarsToFit = (1:8);
        ModelGR_b.fittingOptions.modelVarsToFit = (1:5);

        ModelGR_a = ModelGR_a.formPropensitiesGeneral('ODE_GR_a');
        ModelGR_b = ModelGR_b.formPropensitiesGeneral('ODE_GR_b');
         
        
        %% Load data into the model
        % Load data for GR-alpha
        ModelGR_a_ode = ModelGR_a.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR_a','normgrnuc';'cytGR_a','normgrcyt'},{'Dex_Conc','100'});

        % Load data for GR-beta
        ModelGR_b_ode = ModelGR_b.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR_b','normgrnuc';'cytGR_b','normgrcyt'});

        %% Switch solver to ODE and generate model codes
        ModelGR_a_ode.solutionScheme = 'ode';  % TODO: loop (right now just solve for 100nM Dex)
        ModelGR_b_ode.solutionScheme = 'ode';
        ModelGR_a_ode.useHybrid = false;
        ModelGR_b_ode.useHybrid = false;
        
        %% Solve and make plots
        ODEsoln_a = ModelGR_a_ode.solve; % %TODO: fix broken
        ODEsoln_b = ModelGR_b_ode.solve;
        %plotODEresults(ODEsoln_a,ODEsoln_a,ModelGRfit_a{3})  % These plots only work if
        %plotODEresults(ODEsoln_b,ODEsoln_b,ModelGRfit_b{1})  % FSP case was already run.
                                                              %% Also "plotODEresults" currently 
                                                              %% only works for "extendedModel"!

        %% Fit new parameters to match all ODE data.
        ModelGR_a_ode.fittingOptions.modelVarsToFit = 1:8;
        ModelGR_a_ode.fittingOptions.logPrior = [];
        GRpars_a_ode = [ModelGR_a_ode.parameters{ModelGR_a_ode.fittingOptions.modelVarsToFit,2}];
        GRpars_a_ode = ModelGR_a_ode.maximizeLikelihood(GRpars_a_ode,fitOptions);

        ModelGR_b_ode.fittingOptions.modelVarsToFit = 1:5;
        ModelGR_b_ode.fittingOptions.logPrior = [];
        GRpars_b_ode = [ModelGR_b_ode.parameters{ModelGR_b_ode.fittingOptions.modelVarsToFit,2}];
        GRpars_b_ode = ModelGR_b_ode.maximizeLikelihood(GRpars_b_ode,fitOptions);
            
        %% Plot ODE fit results.
        %ModelGR_b_ode.parameters(ModelGR_b_ode.fittingOptions.modelVarsToFit,2) = num2cell(GRpars_b_ode);
        %ODEsoln_b_dataFit = ModelGR_b_ode.solve;
            %plotODEresults(ModelGR_b_ode,ODEsoln_b_dataFit,ModelGRfit_b{1})
            
            %%    STEP 3.F.1. -- Create New Objective Function Combining all of the previous ones.
            % Remove all priors from individual models.
            ModelGRDusp100nM_ext_red.fittingOptions.logPrior = [];
            for i=1:3
                ModelGRfit{1}.fittingOptions.logPrior = [];
                ModelGRfit{1}.tSpan = ModelGRfit{1}.dataSet.times;
            end
            extendedMod.fittingOptions.logPrior = [];
            fullPars = [extendedMod.parameters{[1:12,14,15],2}];
            % otherewise, we will use the set that was saved in the data dump from
            % the older workspace.

        %%    STEP 3.F.2. -- Fit all objective functions at once.
        % Create prior for all parameters
        log10PriorMean = [-1 -1 0 -2,... %dusp1 pars
            -1 -3 -2 -1 -2 -2 -2 0.5, ...%GR pars
            NaN, ... % Dex concentration -- known
            -2, -3]; % dusp1 transport, cyt RNA degradation
        log10PriorStd = 2*ones(1,15);
        logPriorAll = @(x)-sum((log10(x)-log10PriorMean([1:12,14,15])).^2./(2*log10PriorStd([1:12,14,15]).^2));
        % extendedMod.fittingOptions.modelVarsToFit = [1:12,14,15];
        Organization = {ModelGRfit{1},[5:12],[5:12],'computeLikelihood',1;...
            ModelGRfit{2},[5:12],[5:12],'computeLikelihood',1;...
            ModelGRfit{3},[5:12],[5:12],'computeLikelihood',1;...
            ModelGRDusp100nM_ext_red,[1:12,14],[1:13],'computeLikelihood',1;...
            extendedMod,[1:12,14:15],[1:14],'computeLikelihoodODE',0.01};
        Organization = getTotalFitErr(Organization,fullPars,true);
        getTotalFitErr(Organization,fullPars,false)
        objAll = @(x)-getTotalFitErr(Organization,exp(x),false)-logPriorAll(exp(x));
        for jj = 1:fitIters
            fullParsLog = log(fullPars);
            fullPars = exp(fminsearch(objAll,fullParsLog,fitOptions));
        end
        
        %%    STEP 3.F.3. -- Plot all results
        % Plot GR Distribution
        makeGRPlots(combinedGRModel,fullPars(5:12))
        
        % Plot DUSP1 100nm FIT and other PREDICTED distributions
        showCases = [1,1,1,1];
        ModelGRDusp100nM_ext_red.fittingOptions.modelVarsToFit = [1:12,14];
        ModelGRDusp100nM_ext_red.parameters([1:12,14],2) = num2cell(fullPars([1:13]));
        makePlotsDUSP1({ModelGRDusp100nM_ext_red},ModelGRDusp100nM_ext_red,fullPars([1:13]),Dusp1FitCases,showCases)
        
        % Plot ODE Cyt FIT Results at 100nM Dex
        extendedMod.parameters([1:12,14,15],2) = num2cell(fullPars);
        soln100 = extendedMod.solve;
        plotODEresults(extendedMod,soln100,ModelGRfit{3},500)
        set(gcf,'Name','ODE Fits -- 100nM Dex')
        
        % Plot ODE Predictions at other DEX concentrations
        % extendedMod0p3 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        %     {'rna','RNA_DUSP1_nuc'; ...
        %     'rCyt','RNA_DUSP1_cyto'},...
        %     {'Dex_Conc','0.3'});
        % extendedMod0p3.parameters(13,:) = {'Dex0',0.3};
        % 
        % extendedMod1 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        %     {'rna','RNA_DUSP1_nuc'; ...
        %     'rCyt','RNA_DUSP1_cyto'},...
        %     {'Dex_Conc','1'});
        % extendedMod1.parameters(13,:) = {'Dex0',1.0};
        % 
        % extendedMod10 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        %     {'rna','RNA_DUSP1_nuc'; ...
        %     'rCyt','RNA_DUSP1_cyto'},...
        %     {'Dex_Conc','10'});
        % extendedMod10.parameters(13,:) = {'Dex0',10};
        % 
        % plotODEresults(extendedMod1,extendedMod1.solve,ModelGRfit{1},501)
        % set(gcf,'Name','ODE Predictions -- 1.0nM Dex')
        % 
        % plotODEresults(extendedMod10,extendedMod10.solve,ModelGRfit{2},502)
        % set(gcf,'Name','ODE Predictions -- 10nM Dex')

    case 3  
        %% Create SSA model 
        %SSAModel_100.initialCondition = [2;0;round(soln100.ode(1,3:6))'];
        %SSAModel_100.initialTime = SSAModel_100.tSpan(1);
        %% SSAs        
        SSA_ModelGR_a = SSIT;
        SSA_ModelGR_b = SSIT;
        SSA_ModelGR_a.species = {'cytGR_a';'nucGR_a'};
        SSA_ModelGR_b.species = {'cytGR_b';'nucGR_b'};
        SSA_ModelGR_a.initialCondition = [5;1];
        SSA_ModelGR_b.initialCondition = [5;1];
        SSA_ModelGR_a.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'dc * cytGR_a'; 'kn2c * nucGR_a'; 'dn * nucGR_a'};
        SSA_ModelGR_b.propensityFunctions = {'kb1 * cytGR_b'; 'bb1'; 'dc * cytGR_b'; 'kn2c * nucGR_b'; 'dn * nucGR_b'};
        SSA_ModelGR_a.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        SSA_ModelGR_b.stoichiometry = [-1,1,-1,1,0;...
                                    1,0,0,-1,-1];
        SSA_ModelGR_a.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'ba1',14e-5;...
                                'MDex',5;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;'Dex0',100});
        SSA_ModelGR_b.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'kb1',0.01;'bb1',14e-5});
        SSA_ModelGR_a.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        disp('The GR-alpha model is: ')
        SSA_ModelGR_a.summarizeModel
        disp('and the GR-beta model is: ')
        SSA_ModelGR_b.summarizeModel

        SSA_ModelGR_a.solutionScheme = 'SSA';
        SSA_ModelGR_b.solutionScheme = 'SSA';

        SSA_ModelGR_a.fittingOptions.modelVarsToFit = (1:8);
        SSA_ModelGR_b.fittingOptions.modelVarsToFit = (1:5);

        SSA_ModelGR_a = SSA_ModelGR_a.formPropensitiesGeneral('SSA_GR_a');
        SSA_ModelGR_b = SSA_ModelGR_b.formPropensitiesGeneral('SSA_GR_b');

        SSA_ModelGR_a.tSpan = [-500,SSA_ModelGR_a.tSpan];
        % A negative initial time is needed to allow model to equilibrate before
        % starting.  This causes long run times.
        SSA_ModelGR_b.tSpan = [-500,SSA_ModelGR_b.tSpan];
        SSA_ModelGR_a.initialTime = SSA_ModelGR_a.tSpan(1); % Set initial time
        SSA_ModelGR_b.initialTime = SSA_ModelGR_b.tSpan(1); % Set initial time

        SSA_ModelGR_a.useHybrid = false;
        SSA_ModelGR_b.useHybrid = false;
        SSA_ModelGR_a.ssaOptions.useParalel = true;
        SSA_ModelGR_b.ssaOptions.useParalel = true;

        %% Run SSA Simulations
        ssaSoln_GR_a = SSA_ModelGR_a.solve;
        ssaSoln_GR_b = SSA_ModelGR_b.solve;
        
        %% Plot SSA Results for Cytoplasmic Distributions (100nM Dex)
        % Plot Fits for 100nM Dex 
        %makeCytDistPlots(ssaSoln_GR_a,GR_a_fspSoln,600,2:10,1:9,6,2)
        %makeCytDistPlots(ssaSoln_GR_b,GR_b_fspSoln,600,2:10,1:9,6,2)
        
        %% Predict Cyt distributions for other 0.3nM Dex 
        % SSAModel_0p3 = SSAModel_100;
        % SSAModel_0p3.parameters(13,:) = {'Dex0',0.3};
        % SSAModel_0p3.tSpan = [-500,extendedMod0p3.tSpan];
        % ssaSoln_0p3 = SSAModel_0p3.solve;
        % 
        % %% Predict Cyt distributions for other 1.0nM Dex 
        % SSAModel_1 = SSAModel_100;
        % SSAModel_1.parameters(13,:) = {'Dex0',1.0};
        % SSAModel_1.tSpan = [-500,extendedMod1.tSpan];
        % ssaSoln_1 = SSAModel_1.solve;
        % 
        % %%    STEP 3.F.3. -- Predict Cyt distributions for other 10nM Dex 
        % SSAModel_10 = SSAModel_100;
        % SSAModel_10.parameters(13,:) = {'Dex0',10};
        % SSAModel_10.tSpan = [-500,extendedMod10.tSpan];
        % ssaSoln_10 = SSAModel_10.solve;
        % 
        % %% Make resulting plots
        % makeCytDistPlots(ssaSoln_0p3,extendedMod0p3,601,[2:8],[1:7],6,2)
        % makeCytDistPlots(ssaSoln_1,extendedMod1,602,[3:8],[1:6],6,2)
        % makeCytDistPlots(ssaSoln_10,extendedMod10,603,[3:8],[1:6],6,2)
        % 
        % %%
        % %% Save Results for Easier Use in subsequent runs.
        % %parsAll_GR_Dusp1_TS = [extendedMod.parameters{:,2}];
        % %parsAll_GR_Dusp1_TS(16) = parsAllandTS(end);
        % varNames = unique({'ModelGR'
        %     'GRfitCases'
        %     'log10PriorMean'
        %     'log10PriorStd'
        %     'GRpars'
        %     'ModelGRparameterMap'
        %     'ModelGRfit'
        %     'boundGuesses'
        %     'ModelGRDusp100nM'
        %     'GRfitCases'
        %     'log10PriorMean'
        %     'log10PriorStd'
        %     'duspLogPrior'
        %     'DUSP1pars'
        %     'Dusp1FitCases'
        %     'ModelGRfit'
        %     'extendedMod'
        %     'ModelGRDusp100nM_ext_red'
        %     'fullPars'
        %     'fimResults'
        %     'MHResultsGR'
        %     'MHResultsDusp1'
        %     'sensSoln'
        %     });
        % 
        % save('workspace_Feb_]2024',varNames{:})
end

%% Extra Function

function makeCytDistPlots(ssaSoln_100,extendedMod,fignum,timeIndsMod,timeIndsDat,speciesIndMod,speciesIndDat)
    arguments
        ssaSoln_100
        extendedMod
        fignum = 1;
        timeIndsMod = [];
        timeIndsDat = [];
        speciesIndMod = 1;
        speciesIndDat = 1;
    end
    figure(fignum); clf;
    if isempty(timeIndsMod)
        timeIndsMod = [1:length(ssaSoln_100.T_array)];
    end
    if isempty(timeIndsDat)
        timeIndsDat = [1:length(extendedMod.dataSet.times)];
    end
    
    if length(timeIndsMod)~=length(timeIndsDat)
        error('Length of data and model time points must be equal')
    end
    
    Nrows = ceil(sqrt(length(timeIndsMod)));
    for i = 1:length(timeIndsMod)
        subplot(Nrows,Nrows,i)
    
        % Add SSA to histogram plot
        M = squeeze(ssaSoln_100.trajs(speciesIndMod,timeIndsMod(i),:));
        H = histogram(M,'Normalization','pdf');
        hold on
    
        % Add data to histogram plot
        dMat = double(extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(timeIndsDat(i),:,:));
        N = sum(dMat,"all");
        if speciesIndDat==1
            PD = [0;sum(dMat,2)/N];
        else
            PD = [0;sum(dMat,1)'/N];
        end
        
        binEdges = round(H.BinEdges);
        nBins = length(binEdges)-1;
        PDbinned = zeros(nBins,1);
        binwidth = binEdges(2)-binEdges(1);
        for j = 1:nBins
            PDbinned(j) = sum(PD(binEdges(j)+1:binEdges(j+1)));
        end
    
        PDbinned = [PDbinned;1-sum(PDbinned)];
        stairs(binEdges,PDbinned/binwidth,'linewidth',2)
    
        set(gca,'FontSize',15,'ylim',[0,0.03])
    
    end
end

