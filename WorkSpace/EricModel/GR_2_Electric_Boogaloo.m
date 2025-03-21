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

fitOptions = optimset('Display','iter','MaxIter',200);
 
%%
GR = input('(1) Convolve GR-alpha + GR-beta with FSP;\n(2) ODES for GR-alpha + GR-beta;\n(3) SSAs\nChoose your destiny: ');

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

        data = readtable('../EricModel/EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv');

        % Filter the data if desired for particular Time_index and/or dex_conc
        filteredData = data(data.time == 150 & data.dex_conc == 100, :);

        % Randomly sample 1000 points from the filtered data
        numSamples = min(1000, height(filteredData)); % Ensure we don’t sample more than available data
        randomIndices = randperm(height(filteredData), numSamples);
        normGRnuc = filteredData.normGRnuc(randomIndices);
        normGRcyt = filteredData.normGRcyt(randomIndices);

        % Plot contourf figures and add scatter points
        for g=1:f
            fspsoln_sptensor_a{g} = double(FSPsoln_a.fsp{g}.p.data);
            fspsoln_sptensor_b{g} = double(FSPsoln_b.fsp{g}.p.data);
    
            figure(g)
    
            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a{g}))
            hold on
            scatter(normGRcyt, normGRnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b{g}))
            hold on
            scatter(normGRcyt, normGRnuc, 'r', 'filled') % Red dots
    
            subplot(1,3,3)
            contourf(log10(conv2solnTensor{g}))
            hold on
            scatter(normGRcyt, normGRnuc, 'r', 'filled') % Red dots
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
        ModelGRparameterMap_b = cell(1,size(GRfitCases,1));
        ModelGRfit_a = cell(1,size(GRfitCases,1));
        ModelGRfit_b = cell(1,size(GRfitCases,1));

        % ModelGRODEfit = cell(1,size(GRfitCases,1));
        
        %% Load data, fit 3 Dex concentrations to GR-alpha's 'Dex0' parameter and set solution scheme to FSP.
        for i=1:3
            ModelGRfit_a{i} = ModelGR_a.loadData("../EricModel/EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv",...
                                                {'nucGR_a','normGRnuc';'cytGR_a','normGRcyt'},...
                                                {'dex_conc',GRfitCases{i,2}});
            ModelGRfit_b{i} = ModelGR_b.loadData("../EricModel/EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv",...
                                                {'nucGR_b','normGRnuc';'cytGR_b','normGRcyt'},...
                                                {'dex_conc',GRfitCases{i,2}});
            ModelGRfit_a{i}.parameters(9,:) = {'Dex0', str2num(GRfitCases{i,1})};
            ModelGRparameterMap_a(i) = {(1:8)};  
            ModelGRparameterMap_b(i) = {(1:5)};
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
    fitMHiters = 1; % Adjust as needed
    
    for GR = 1:fitMHiters
        % Specify dataset time points for GR-alpha.    
        for i = 1:3
            ModelGRfit_a{i}.tSpan = ModelGRfit_a{i}.dataSet.times;
            ModelGRfit_b{i}.tSpan = ModelGRfit_b{i}.dataSet.times;
        end

        % Specify dataset time points for GR-beta.
        %ModelGRfit_b{1}.tSpan = ModelGRfit_b{1}.dataSet.times;

        
        %%
        % PA = rand(2); PA = PA/(sum(PA,"all"));
        % PB = rand(2); PB = PB/(sum(PB,"all"));        
        % % PA=[.5 .15;.25 .1];PB =PA;
        % PC = conv2(PA,PB)        
        % %%       
        % N = 100000;
        % PAs = zeros(size(PA));        
        % PBs = zeros(size(PB));        
        % PCs = zeros(size(PA,1)+size(PB,1)-1,...        
        %     size(PA,2)+size(PB,2)-1);        
        % for i = 1:N        
        %     r = rand;        
        %     j = 1;        
        %     while sum(PA(1:j))<r
        %         j = j+1;
        %     end
        %     PAs(j) = PAs(j) + 1/N;
        %     [kA1,kA2] = ind2sub(size(PAs),j);        
        %     r = rand;
        %     j = 1;
        %     while sum(PB(1:j))<r        
        %         j = j+1;       
        %     end        
        %     PBs(j) = PBs(j) + 1/N;
        %     [kB1,kB2] = ind2sub(size(PBs),j);
        %     PCs(kA1+kB1-1,kA2+kB2-1) = PCs(kA1+kB1-1,kA2+kB2-1)+1/N;
        % end       
        % PCs        
        % PC
        % 
        % %% REMOVE
        % ModelGR_a.solutionScheme = 'FSP';
        % FSPsoln = ModelGR_a.solve;             
        % solnTensor = double(FSPsoln.fsp{3}.p.data);
        % contourf(log10(solnTensor));
        % 
        % %% sample fake data        
        % N = 100;        
        % PAs = zeros(size(solnTensor));        
        % for i = 1:N      
        %     r = rand;        
        %     j = 1;       
        %     while sum(solnTensor(1:j))<r        
        %         j = j+1;       
        %     end
        %     PAs(j) = PAs(j)+1;        
        %     % [kA1,kA2] = ind2sub(size(PAs),j);        
        % end
        % 
        % %% add scatter points to plot.
        % hold on
        % [rows,cols,vals] = find(PAs);        
        % scatter(rows,cols,20,vals,'r','filled')
        % 
        % %% compute likelihood
        % logLikelihood = sum(log(solnTensor(PAs>0)).*PAs((PAs>0)),"all");


        %% Solve for combined GR-alpha model for three dex_conc=Dex0 and fit using a single parameter set.
        % for jj = 1:fitIters
        %     % Solve
        %     combinedGRModel_a = SSITMultiModel(ModelGRfit_a,ModelGRparameterMap_a,logPriorGR_a);
        %     %combinedGRModel_a = combinedGRModel_a.initializeStateSpaces(boundGuesses);
        %     [combinedGRModel_a,combinedGRModel_b,conv2solnTensor] = initializeStateSpaces(combinedGRModel_a,combinedGRModel_b,boundGuesses);
        %     %combinedGRModel_a = combinedGRModel_a.updateModels(GRpars_a,false);
        %     %GRpars_a = combinedGRModel_a.maximizeLikelihood(GRpars_a, fitOptions);
        %     save('combinedGRModel_a','GRpars_a') 
        % end
        
        %% Solve for combined GR-alpha model for three dex_conc=Dex0 and fit using a single parameter set.
        for jj = 1:fitIters
            % Solve
            combinedGRModel_a = SSITMultiModel(ModelGRfit_a,ModelGRparameterMap_a,logPriorGR_a);
            combinedGRModel_b = SSITMultiModel(ModelGRfit_b,ModelGRparameterMap_b,logPriorGR_b);
            [combinedGRModel_a,combinedGRModel_b,conv2solnTensor] = solveandconvolve(combinedGRModel_a,combinedGRModel_b,boundGuesses);
            combinedGRModel_a = combinedGRModel_a.updateModels(GRpars_a,false);
            combinedGRModel_b = combinedGRModel_b.updateModels(GRpars_b,false);
            %combinedGRModel_a = combinedGRModel_a.initializeStateSpaces(boundGuesses);
            %combinedGRModel_a = combinedGRModel_a.updateModels(GRpars_a,false);
            %GRpars_a = combinedGRModel_a.maximizeLikelihood(GRpars_a, fitOptions);
            save('combinedGRModel_a','GRpars_a') 
            %% Solve for GR-beta model.
            combinedGRModel_b = SSITMultiModel(ModelGRfit_b,ModelGRparameterMap_b,logPriorGR_b);
            %combinedGRModel_b = combinedGRModel_b.initializeStateSpaces(boundGuesses);
            %combinedGRModel_b = combinedGRModel_b.updateModels(GRpars_b,false);
            %GRpars_b = combinedGRModel_b.maximizeLikelihood(GRpars_b, fitOptions);
            save('combinedGRModel_b','GRpars_b') 

            [combinedGRModel_a,combinedGRModel_b,conv2solnTensor] = solveandconvolve(combinedGRModel_a,combinedGRModel_b,boundGuesses);

        end
        
        %% Solve for GR-beta model.
        % ModelGRfit_b{1}.fspOptions.fspTol = 1e-4;
        % ModelGRfit_b{1}.fspOptions.bounds = boundGuesses{1};
        % [GR_b_fspSoln,ModelGRfit_b{1}.fspOptions.bounds] = ModelGRfit_b{1}.solve;
        % [GR_b_fspSoln,ModelGRfit_b{1}.fspOptions.bounds] = ModelGRfit_b{1}.solve(GR_b_fspSoln.stateSpace);
        % GRpars_b = ModelGRfit_b{1}.maximizeLikelihood(GRpars_b, fitOptions);
        % save('ModelGRfit_b','GRpars_b','GR_b_fspSoln')             
        % 
        % save('GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln')
        

        %% Compute FIM 
        % for ModelGR_a parameters.
        combinedGRModel_a = combinedGRModel_a.computeFIMs([],'log');

        fimGR_a_withPrior = combinedGRModel_a.FIM.totalFIM+... % the FIM in log space.
            diag(1./(log10PriorStd_a(ModelGR_a.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.

        % for ModelGR_b parameters.
        combinedGRModel_b = combinedGRModel_b.computeFIMs([],'log');

        fimGR_b_withPrior = combinedGRModel_b.FIM.totalFIM+... % the FIM in log space.
            diag(1./(log10PriorStd_b(ModelGR_b.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.
        %ModelGR_b_fimResults = ModelGRfit_b{1}.computeFIM([],'log'); % Compute individual FIMs

        %fimGR_b_withPrior = ModelGRfit_b{1}.totalFim(ModelGR_b_fimResults,ModelGRfit_b{1}.dataSet.nCells,diag(1./(log10PriorStd_b(ModelGRfit_b{1}.fittingOptions.modelVarsToFit)*log(10)).^2));  % Add prior in log space.

        %
        if min(eig(fimGR_a_withPrior))<1
            disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
            fimGR_a_withPrior = fimGR_a_withPrior + 1*eye(length(fimGR_a_withPrior));
        end
        fimGR_a_covFree = fimGR_a_withPrior^-1;
        fimGR_a_covFree = 0.5*(fimGR_a_covFree+fimGR_a_covFree');

        if min(eig(fimGR_b_withPrior))<1
             disp('Warning -- FIM has one or more small eigenvalues. Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
             fimGR_b_withPrior = fimGR_b_withPrior + 1*eye(length(fimGR_b_withPrior));
        end
        fimGR_b_covFree = fimGR_b_withPrior^-1;
        fimGR_b_covFree = 0.5*(fimGR_b_covFree+fimGR_b_covFree');
    
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
        % MHfitOptions_b.thin=1;
        % MHfitOptions_b.numberOfSamples=1000;
        % MHfitOptions_b.burnIn=100;
        % MHfitOptions_b.progress=true;
        % MHfitOptions_b.numChains = 1;
        % 
        % % Use FIM computed above rather than making SSIT call 'useFIMforMetHast'
        % % which forces SSIT.m to compute it within (no prior, etc.)
        % MHfitOptions_b.useFIMforMetHast = false;
        % MHfitOptions_b.proposalDistribution=@(x)mvnrnd(x,fimGR_b_covFree);
        % 
        % MHfitOptions_b.saveFile = 'TMPEricMHGR_b.mat';
        % [~,~,MHResultsGR_b] = combinedGRModel_b.maximizeLikelihood(...
        %     GRpars_b, MHfitOptions_b, 'MetropolisHastings');
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
    data = readtable('../EricModel/EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv'); 

    tspan = [0,10,30,50,75,150,180];

    for model=1:size(GRfitCases,1)
        dataDex = data(data.dex_conc == str2double(GRfitCases{model}), :);
        %for t=1:max(numel(fspSolnsSMM(model).fsp),numel(GR_b_fspSoln.fsp))
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
            normGRnuc = dataDexTime.normGRnuc(randomIndices);
            normGRcyt = dataDexTime.normGRcyt(randomIndices);

            % t is time point, so solution tensors are FSP probabilities across states for each time point
            conv2solnTensor_postData{t} = conv2(double(fspSolnsSMM(model).fsp{t}.p.data),double(GR_b_fspSoln.fsp{t}.p.data));
            figure(figIndex)
            contourf(log10(conv2solnTensor_postData{t}))
            hold on 
            scatter(normGRcyt, normGRnuc, 'r', 'filled')

            % Add subplots (can be commented out)
            fspsoln_sptensor_a_postData{t} = double(fspSolnsSMM(model).fsp{t}.p.data);
            fspsoln_sptensor_b_postData{t} = double(GR_b_fspSoln.fsp{t}.p.data);

            subplot(1,3,1)
            contourf(log10(fspsoln_sptensor_a_postData{t}))
            hold on
            scatter(normGRcyt, normGRnuc, 'r', 'filled') % Red dots

            subplot(1,3,2)
            contourf(log10(fspsoln_sptensor_b_postData{t}))
            hold on
            scatter(normGRcyt, normGRnuc, 'r', 'filled') % Red dots

            subplot(1,3,3)
            contourf(log10(conv2solnTensor_postData{t}))
            hold on
            scatter(normGRcyt, normGRnuc, 'r', 'filled') % Red dots
        end
    end

    %% resample data if wanted

    data = readtable('../EricModel/EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv');

        % Filter the data if desired for particular Time_index and/or dex_conc
        filteredData = data(data.time == 180 & data.dex_conc == 10, :);

        % Randomly sample 1000 points from the filtered data
        numSamples = min(1000, height(filteredData)); % Ensure we don’t sample more than available data
        randomIndices = randperm(height(filteredData), numSamples);
        normGRnuc = filteredData.normGRnuc(randomIndices);
        normGRcyt = filteredData.normGRcyt(randomIndices);

    %% compute likelihood
    
        normgr = [filteredData.normGRcyt,filteredData.normGRnuc];
        for logL=1:max(size(conv2solnTensor_postData))
            logLikelihood = sum(log(conv2solnTensor_postData{logL}(normgr>0)).*normgr((normgr>0)),"all")
        end
    %%     STEP 1.F. -- Make Plots of GR Fit Results
        makeGRPlots(combinedGRModel_a,GRpars_a)

        % ModelGRfit_b{1}.makeFitPlot([],1,[],true,'STD')

        makeGRPlots(combinedGRModel_b,GRpars_b)
    

    save('conv2solnTensor_postData','GRpars_a','combinedGRModel_a', 'ModelGRfit_a','ModelGRfit_b','GRpars_b','GR_b_fspSoln', 'fspSolnsSMM')

    case 2 
        %% ODEs        
        ODE_ModelGR = SSIT;
        ODE_ModelGR.species = {'cytGR_a';'nucGR_a';'cytGR_b';'nucGR_b'};
        ODE_ModelGR.initialCondition = [7;2;8;2];
        ODE_ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'dc * cytGR_a'; 'kn2c * nucGR_a'; 'dn * nucGR_a';... 
                                            'kcn * cytGR_b'; 'bb1'; 'dc * cytGR_b'; 'kn2c * nucGR_b'; 'dn * nucGR_b'};
        ODE_ModelGR.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0;...
                                      1,0,0,-1,-1,0,0,0,0,0;...
                                      0,0,0,0,0,-1,1,-1,1,0;...
                                      0,0,0,0,0,1,0,0,-1,-1];
        ODE_ModelGR.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'ba1',14e-5;...
                                   'MDex',5;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;...
                                   'kcn',0.01;'bb1',14e-5;'Dex0',100});
        ODE_ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};

        ODE_ModelGR.solutionScheme = 'ode';
        ODE_ModelGR.useHybrid = false;
        ODE_ModelGR.fittingOptions.modelVarsToFit = (1:10);

        disp('The GR model is: ')
        ODE_ModelGR.summarizeModel

        ODE_ModelGR = ODE_ModelGR.formPropensitiesGeneral('ODE_GR');

        %% Solve ODE and make plots
        ODE_GR_soln = ODE_ModelGR.solve; 
        plotODE(ODE_GR_soln,ODE_ModelGR.species)

        dusp1 = input('(1) GR-beta turns off the DUSP1 gene;\n(2) GR-beta has no effect on DUSP1;\nChoose your destiny: ');

        switch dusp1
            case 1 
                %% GR-beta turns off DUSP1
                % Add DUSP1 to model
                ODE_ModelDUSP = ODE_ModelGR;
                ODE_ModelDUSP = ODE_ModelDUSP.addSpecies({'offGene'},2);
                ODE_ModelDUSP = ODE_ModelDUSP.addSpecies({'onGene'},0);
                ODE_ModelDUSP = ODE_ModelDUSP.addSpecies({'rna'},5);
        
                ODE_ModelDUSP.propensityFunctions{11,1} = 'kon*offGene*nucGR_a';
                ODE_ModelDUSP.propensityFunctions{12,1} = 'koff*onGene*nucGR_b';
                ODE_ModelDUSP.propensityFunctions{13,1} = 'kr*onGene';
                ODE_ModelDUSP.propensityFunctions{14,1} = 'dr*rna';
        
                ODE_ModelDUSP.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0,0,0,0,0;...
                                                1,0,0,-1,-1,0,0,0,0,0,0,0,0,0;...
                                                0,0,0,0,0,-1,1,-1,1,0,0,0,0,0;...
                                                0,0,0,0,0,1,0,0,-1,-1,0,0,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,-1,1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,1,-1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,0,0,1,-1];
        
                ODE_ModelDUSP.parameters(12,:) = {'koff',0.1};
                ODE_ModelDUSP.parameters(13,:) = {'kon',0.1};
                ODE_ModelDUSP.parameters(14,:) = {'kr',1};
                ODE_ModelDUSP.parameters(15,:) = {'dr',0.02};
        
                ODE_ModelDUSP.summarizeModel
        
                %% Solve ODE
                ODE_DUSP_soln = ODE_ModelDUSP.solve;
        
                %% Plot ODE Results (100nM Dex) 
                plotODE(ODE_DUSP_soln,ODE_ModelDUSP.species)
            case 2
                %% GR-beta has no effect on DUSP1
                % Add DUSP1 to model
                ODE_ModelDUSP = ODE_ModelGR;
                ODE_ModelDUSP = ODE_ModelDUSP.addSpecies({'offGene'},2);
                ODE_ModelDUSP = ODE_ModelDUSP.addSpecies({'onGene'},0);
                ODE_ModelDUSP = ODE_ModelDUSP.addSpecies({'rna'},5);
        
                ODE_ModelDUSP.propensityFunctions{11,1} = 'kon*offGene*nucGR_a';
                ODE_ModelDUSP.propensityFunctions{12,1} = 'koff*onGene';
                ODE_ModelDUSP.propensityFunctions{13,1} = 'kr*onGene';
                ODE_ModelDUSP.propensityFunctions{14,1} = 'dr*rna';
        
                ODE_ModelDUSP.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0,0,0,0,0;...
                                                1,0,0,-1,-1,0,0,0,0,0,0,0,0,0;...
                                                0,0,0,0,0,-1,1,-1,1,0,0,0,0,0;...
                                                0,0,0,0,0,1,0,0,-1,-1,0,0,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,-1,1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,1,-1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,0,0,1,-1];
        
                ODE_ModelDUSP.parameters(12,:) = {'koff',0.1};
                ODE_ModelDUSP.parameters(13,:) = {'kon',0.1};
                ODE_ModelDUSP.parameters(14,:) = {'kr',1};
                ODE_ModelDUSP.parameters(15,:) = {'dr',0.02};
        
                ODE_ModelDUSP.summarizeModel
        
                %% Solve ODE
                ODE_DUSP_soln = ODE_ModelDUSP.solve;
        
                %% Plot ODE Results (100nM Dex) 
                plotODE(ODE_DUSP_soln,ODE_ModelDUSP.species)
        end
        
        %% Load data into the model
        % Load data for GR-alpha
        ModelGR_ode = ModelGR.loadData("../EricModel/EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv",...
                {'nucGR_a','normGRnuc';'cytGR_a','normGRcyt'},{'Dex_Conc','100'});
        
        %% Solve and make plots
        ODEsoln_GR = ModelGR_ode.solve; % %TODO: fix broken
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
        %     {'dex_conc','0.3'});
        % extendedMod0p3.parameters(13,:) = {'Dex0',0.3};
        % 
        % extendedMod1 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        %     {'rna','RNA_DUSP1_nuc'; ...
        %     'rCyt','RNA_DUSP1_cyto'},...
        %     {'dex_conc','1'});
        % extendedMod1.parameters(13,:) = {'Dex0',1.0};
        % 
        % extendedMod10 = extendedMod.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        %     {'rna','RNA_DUSP1_nuc'; ...
        %     'rCyt','RNA_DUSP1_cyto'},...
        %     {'dex_conc','10'});
        % extendedMod10.parameters(13,:) = {'Dex0',10};
        % 
        % plotODEresults(extendedMod1,extendedMod1.solve,ModelGRfit{1},501)
        % set(gcf,'Name','ODE Predictions -- 1.0nM Dex')
        % 
        % plotODEresults(extendedMod10,extendedMod10.solve,ModelGRfit{2},502)
        % set(gcf,'Name','ODE Predictions -- 10nM Dex')

    case 3  
        %% Create SSA model of GR
        %SSAModel_100.initialCondition = [2;0;round(soln100.ode(1,3:6))'];
        %SSAModel_100.initialTime = SSAModel_100.tSpan(1);
        %% SSAs        
        SSA_ModelGR = SSIT;
        SSA_ModelGR.species = {'cytGR_a';'nucGR_a';...
                               'cytGR_b';'nucGR_b'};
        SSA_ModelGR.initialCondition = [7;2;8;2];
        SSA_ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'dc * cytGR_a'; 'kn2c * nucGR_a'; 'dn * nucGR_a';... 
                                            'kcn * cytGR_b'; 'bb1'; 'dc * cytGR_b'; 'kn2c * nucGR_b'; 'dn * nucGR_b'};
        SSA_ModelGR.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0;...
                                      1,0,0,-1,-1,0,0,0,0,0;...
                                      0,0,0,0,0,-1,1,-1,1,0;...
                                      0,0,0,0,0,1,0,0,-1,-1];
        SSA_ModelGR.parameters = ({'kn2c',0.01;'dc',1e-5;'dn',1e-6;'ba1',14e-5;...
                                   'MDex',5;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;...
                                   'kcn',0.01;'bb1',14e-5;'Dex0',100});
        SSA_ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        disp('The GR model is: ')
        SSA_ModelGR.summarizeModel

        SSA_ModelGR.solutionScheme = 'SSA';

        SSA_ModelGR.fittingOptions.modelVarsToFit = (1:10);

        SSA_ModelGR = SSA_ModelGR.formPropensitiesGeneral('SSA_GR');

        SSA_ModelGR.tSpan = [-500,SSA_ModelGR.tSpan];
        % A negative initial time is needed to allow model to equilibrate before
        % starting.  This causes long run times.
        SSA_ModelGR.initialTime = SSA_ModelGR.tSpan(1); % Set initial time

        SSA_ModelGR.useHybrid = false;
        SSA_ModelGR.ssaOptions.useParalel = true;

        %% Run SSA Simulations
        ssaSoln_GR = SSA_ModelGR.solve;
        
        %% Plot SSA Results (100nM Dex) 
        plotSSA(ssaSoln_GR, 'all', 2200);

        dusp1 = input('(1) GR-beta turns off the DUSP1 gene;\n(2) GR-beta has no effect on DUSP1;\nChoose your destiny: ');

        switch dusp1
            case 1 
                %% GR-beta turns off DUSP1
                % Add DUSP1 to model
                SSA_ModelDUSP = SSA_ModelGR;
                SSA_ModelDUSP = SSA_ModelDUSP.addSpecies({'offGene'},2);
                SSA_ModelDUSP = SSA_ModelDUSP.addSpecies({'onGene'},0);
                SSA_ModelDUSP = SSA_ModelDUSP.addSpecies({'rna'},5);
        
                SSA_ModelDUSP.propensityFunctions{11,1} = 'kon*offGene*nucGR_a';
                SSA_ModelDUSP.propensityFunctions{12,1} = 'koff*onGene*nucGR_b';
                SSA_ModelDUSP.propensityFunctions{13,1} = 'kr*onGene';
                SSA_ModelDUSP.propensityFunctions{14,1} = 'dr*rna';
        
                SSA_ModelDUSP.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0,0,0,0,0;...
                                                1,0,0,-1,-1,0,0,0,0,0,0,0,0,0;...
                                                0,0,0,0,0,-1,1,-1,1,0,0,0,0,0;...
                                                0,0,0,0,0,1,0,0,-1,-1,0,0,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,-1,1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,1,-1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,0,0,1,-1];
        
                SSA_ModelDUSP.parameters(12,:) = {'koff',0.1};
                SSA_ModelDUSP.parameters(13,:) = {'kon',0.1};
                SSA_ModelDUSP.parameters(14,:) = {'kr',1};
                SSA_ModelDUSP.parameters(15,:) = {'dr',0.02};
        
                SSA_ModelDUSP.summarizeModel
        
                %% Run SSA Simulations
                ssaSoln_DUSP = SSA_ModelDUSP.solve;
        
                %% Plot SSA Results (100nM Dex) 
                plotSSA(ssaSoln_DUSP, 'all', 2200);
            case 2
                %% GR-beta has no effect on DUSP1
                % Add DUSP1 to model
                SSA_ModelDUSP = SSA_ModelGR;
                SSA_ModelDUSP = SSA_ModelDUSP.addSpecies({'offGene'},2);
                SSA_ModelDUSP = SSA_ModelDUSP.addSpecies({'onGene'},0);
                SSA_ModelDUSP = SSA_ModelDUSP.addSpecies({'rna'},5);
        
                SSA_ModelDUSP.propensityFunctions{11,1} = 'kon*offGene*nucGR_a';
                SSA_ModelDUSP.propensityFunctions{12,1} = 'koff*onGene';
                SSA_ModelDUSP.propensityFunctions{13,1} = 'kr*onGene';
                SSA_ModelDUSP.propensityFunctions{14,1} = 'dr*rna';
        
                SSA_ModelDUSP.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0,0,0,0,0;...
                                                1,0,0,-1,-1,0,0,0,0,0,0,0,0,0;...
                                                0,0,0,0,0,-1,1,-1,1,0,0,0,0,0;...
                                                0,0,0,0,0,1,0,0,-1,-1,0,0,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,-1,1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,1,-1,0,0;...
                                                0,0,0,0,0,0,0,0,0,0,0,0,1,-1];
        
                SSA_ModelDUSP.parameters(12,:) = {'koff',0.1};
                SSA_ModelDUSP.parameters(13,:) = {'kon',0.1};
                SSA_ModelDUSP.parameters(14,:) = {'kr',1};
                SSA_ModelDUSP.parameters(15,:) = {'dr',0.02};
        
                SSA_ModelDUSP.summarizeModel
        
                %% Run SSA Simulations
                ssaSoln_DUSP = SSA_ModelDUSP.solve;
        
                %% Plot SSA Results (100nM Dex)  
                plotSSA(ssaSoln_DUSP, 'all', 2200, SSA_ModelDUSP.species);
        end
end

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
                 {'dex_conc','100'});
             ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
             ModelDusp1parameterMap{i} = (1:4);
             % Set Dex concentration.
             ModelDusp1Fit{i}.parameters{9,2} = str2num(Dusp1FitCases{i,1});
             ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
        end
        DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];

    %% Load data
        ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
                                                        {'rna','RNA_DUSP1_nuc'},{'dex_conc','100'});
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
                 {'dex_conc','100'});
             ModelDusp1Fit{i}.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
             ModelDusp1parameterMap{i} = (1:4);
             % Set Dex concentration.
             ModelDusp1Fit{i}.parameters{13,2} = str2num(Dusp1FitCases{i,1});
             ModelDusp1Fit{i} = ModelDusp1Fit{i}.formPropensitiesGeneral(['EricModDusp1_',num2str(i),'_FSP']);
        end
        DUSP1pars = [ModelDusp1Fit{i}.parameters{ModelGRDusp100nM.fittingOptions.modelVarsToFit,2}];

    %%
    ModelGRDusp100nM = ModelGRDusp100nM.loadData('EricData/pdoCalibrationData_EricIntensity_DexSweeps.csv',...
        {'rna','totalNucRNA'},{'dex_conc','100'});
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
%logLikelihood = sum(log(conv2solnTensor_postData{t}(PAs>0)).*PAs((PAs>0)),"all");

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

