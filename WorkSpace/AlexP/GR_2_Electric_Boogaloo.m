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
    
 
%%
GR = input('(0) FSP for GR-alpha only (base model);\n(1) Convolve GR-alpha + GR-beta;\n(2) ODES for GR-alpha + GR-beta.\nChoose your destiny: ');
%GR = input('(0) GR-alpha only (base model);\n(0) PDO (GR-beta treated as distortion to GR-alpha data);\n(1) Convolve P(GR-alpha) and P(GR-beta);\n(1) The whole kit & caboodle (GR-alpha + GR-beta).\nChoose your destiny: ');

switch GR
    case 0
    % GR-alpha setup    
        disp('You have chosen the base model, GR-alpha.')
        %% STEP 0.B. -- Create Base Model for GR Only
        %    Here, I set up the model for the GR translocation dynamics.
        fitOptions = optimset('Display','iter','MaxIter',300);
    
        % Create blank SSIT model.
        ModelGR = SSIT;
    
        % Set species names.
        ModelGR.species = {'cytGR';'nucGR'};
    
        % Set initial condition.
        ModelGR.initialCondition = [20;1];
    
        % Define propensity functions and input signals:
        ModelGR.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';'knc*nucGR';...
            'kg1';'gg1*cytGR';'gg2*nucGR'};
    
        % Define stoichiometry.
        ModelGR.stoichiometry = [-1,1,1,-1,0;...
            1,-1,0,0,-1];
    
        % Specify parameter guesses.
        ModelGR.parameters = ({'koff',0.1;'kon',0.1;'kr',1;'gr',0.02;...
            'kcn0',0.005;'kcn1',0.02;'gDex',0.003;'knc',0.01;'kg1',14e-5;...
            'gg1',1e-5;'gg2',1e-6;'MDex',5;'Dex0',100});
    
        % Print visual summary of the model
        ModelGR.summarizeModel
    
        %% The log prior will be applied to the fit to multiple models as an additional constraint.
        log10PriorMean = [-1 -1 0 -2,...
            -1 -3 -2 -1 -2 -2 -2 0.5 2];
        log10PriorStd = 2*ones(1,13);
        
        % So it is left out of the prior, since we only want it to be calculated once.
        ModelGR.fittingOptions.logPrior = [];  
    
        
        ModelGR.fspOptions.initApproxSS = true;
        ModelGR.fittingOptions.modelVarsToFit = (5:12);
        
        ModelGR.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        ModelGR = ModelGR.formPropensitiesGeneral('EricRonModGR');
        ModelGR.customConstraintFuns = {'cytGR+nucGR'};
        % TODO - Alex - Constraints for removing stiff dimensions.
        
        [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve;
        [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve(FSPGrSoln.stateSpace);
    
        % STEP 0.B.2. -- Define GR parameters
        GRpars = cell2mat(ModelGR.parameters(5:12,2))';  
    
        % STEP 0.B.3. -- Associate GR Data with Different Instances of Model (10,100nm Dex)
        GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
            '10','10',102,'GR Fit (10nM Dex)';...
            '100','100',103,'GR Fit (100nM Dex)'};
        ModelGRparameterMap = cell(1,size(GRfitCases,1));
        ModelGRfit = cell(1,size(GRfitCases,1));
        % ModelGRODEfit = cell(1,size(GRfitCases,1));
        for i=1:3
            ModelGRfit{i} = ModelGR.loadData("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                {'nucGR','normgrnuc';'cytGR','normgrcyt'},...
                {'Dex_Conc',GRfitCases{i,2}});
            ModelGRfit{i}.parameters(13,:) = {'Dex0', str2num(GRfitCases{i,1})};
            ModelGRparameterMap(i) = {(1:8)};
            % parameters 1 - 8 refer to the parameter set that is relevant to
            % the entire class of models.  In this case, these refer to
            % global parameters 5:13 (GR parameters).
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
    fitIters = 3;
    fitMHiters = 2;
    
    for GR = 1:fitMHiters
        % STEP 1.A. -- Specify dataset time points.    
        for i = 1:3
            ModelGRfit{i}.tSpan = ModelGRfit{i}.dataSet.times;
        end
    
        % STEP 1.B. -- Specify log prior (NOTE: must transpose due to Matlab update that
        %     no longer correctly assumes format when adding single value vector to
        %     column vector).
    
        logPriorGR = @(x)-sum((log10(x)-log10PriorMean(5:12)').^2./(2*log10PriorStd(5:12)'.^2));
    
        % STEP 1.C. -- Combine all three GR models and fit using a single parameter set.
        for jj = 1:fitIters
            combinedGRModel = SSITMultiModel(ModelGRfit,ModelGRparameterMap,logPriorGR);
            combinedGRModel = combinedGRModel.initializeStateSpaces(boundGuesses);
            combinedGRModel = combinedGRModel.updateModels(GRpars,false);
            GRpars = combinedGRModel.maximizeLikelihood(...
                GRpars, fitOptions);
            save('EricModel_MMDex','GRpars') 
        end
    
        save('EricModelGR_MMDex','GRpars','combinedGRModel', 'ModelGRfit', 'log10PriorStd') 
    
        %% STEP 1.D. -- Compute FIM for GR parameters.
        combinedGRModel = combinedGRModel.computeFIMs([],'log');
        fimGR_withPrior = combinedGRModel.FIM.totalFIM+... % the FIM in log space.
            diag(1./(log10PriorStd(ModelGR.fittingOptions.modelVarsToFit)*log(10)).^2);  % Add prior in log space.
    
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
    
    %%     STEP 1.F. -- Make Plots of GR Fit Results
    makeGRPlots(combinedGRModel,GRpars)
    
    save('EricModelGR_MMDex','GRpars','combinedGRModel','MHResultsGR') 
    save('workspaceDec9_2024.mat','GRpars', 'ModelGRfit', 'combinedGRModel','MHResultsGR', 'log10PriorStd')

    case 1
    %% GR-beta setup
        disp('You have chosen GR-alpha + GR-beta by convolution.  WARNING: UNDER CONSTRUCTION')
        ModelGR_b = SSIT;
        ModelGR_b.species = {'cytGR_a';'nucGR_a';'cytGR_b';'nucGR_b'};
        ModelGR_b.initialCondition = [20;1;10;11];
        ModelGR_b.propensityFunctions = {'(kcn0 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR_a'; 'ba1'; 'da1 * cytGR_a'; 'ka1 * nucGR_a'; 'da2 * nucGR_a';...
                                            'bb1'; 'db1 * cytGR_b'; 'kb1 * cytGR_b'; 'kb2 * nucGR_b'; 'db2 * nucGR_b'};
        ModelGR_b.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0;...
                                    1,0,0,-1,-1,0,0,0,0,0;...
                                    0,0,0,0,0,1,-1,-1,1,0;...
                                    0,0,0,0,0,0,0,1,-1,-1];
        ModelGR_b.parameters = ({'MDex',5;'Dex0',100;'gDex',0.003;'kcn0',0.005;'kcn1',0.02;...
                                'ka1',0.01;'ba1',14e-5;'da1',1e-5;'da2',1e-6; ...
                                'kb1',0.01;'kb2',0.01;'bb1',14e-5;'db1',1e-5;'db2',1e-6});
        ModelGR_b.inputExpressions = {'IDex','Dex0*exp(-gDex*t)'};
        ModelGR_b.summarizeModel

        % The log prior will be applied to the fit to multiple models as an additional constraint.
        %log10PriorMean_b = [0.5 2 -3 -2 -2 -2 -4 -5 -6 -1 -1 -4 -5 -6];
        %log10PriorStd_b = 2*ones(1,14);

    case 2
        %% Solve ODEs for fancy GR models
        disp('You have chosen GR-alpha + GR-beta by ODEs.')
        
        % Load experimental data
        data = readtable("../EricModel/EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv");  
        %data = readtable("~/testGR_ode.txt"); % test 
        
        % Extract unique cell IDs
        cell_ids = unique(data.Cell_id);
        
        % Initialize storage for simulation results
        all_simulations = table();
        
        %for i = 1:length(cell_ids)
        for i = 1:100 % WARNING:  hacky sack!
            % Extract data for the current cell
            cell_data = data(data.Cell_id == cell_ids(i), :);
        
            % Extract time and measurements
            time_data = cell_data.time;
            normgrnuc = cell_data.normgrnuc; % nucGR_a + nucGR_b
            normgrcyt = cell_data.normgrcyt; % cytGR_a + cytGR_b
        
            % Initial conditions for this cell
            x1_0 = 20.0; % Initial concentration of cytGR_a
            x2_0 = 1.0; % Initial concentration of nucGR_a
            x3_0 = 10.0; % Initial concentration of cytGR_b
            x4_0 = 1.0; % Initial concentration of nucGR_b
            y0 = [x1_0; x2_0; x3_0; x4_0];
        
            % Time span for this cell
            tspan = [min(time_data), max(time_data)];
        
            % Solve the ODEs
            [t, y] = ode45(@crs_odes, tspan, y0);
        
            % Extract results
            x1 = y(:, 1); % cytGR_a
            x2 = y(:, 2); % nucGR_a
            x3 = y(:, 3); % cytGR_b
            x4 = y(:, 4); % nucGR_b
        
            % Compute aggregate values for comparison
            sim_normgrnuc = x2 + x4; % nucGR_a + nucGR_b
            sim_normgrcyt = x1 + x3; % cytGR_a + cytGR_b
        
            % Store results in a table
            sim_table = table(t, sim_normgrnuc, sim_normgrcyt, repmat(cell_ids(i), length(t), 1), ...
                              'VariableNames', {'time', 'sim_normgrnuc', 'sim_normgrcyt', 'Cell_id'});
            
            % Append to all simulations
            all_simulations = [all_simulations; sim_table];
        
            % Plot results for the current cell
            figure;
            plot(t, sim_normgrnuc, 'LineWidth', 2, 'DisplayName', 'Simulated normgrnuc');
            hold on;
            plot(t, sim_normgrcyt, 'LineWidth', 2, 'DisplayName', 'Simulated normgrcyt');
            plot(time_data, normgrnuc, 'o', 'DisplayName', 'Experimental normgrnuc');
            plot(time_data, normgrcyt, 'o', 'DisplayName', 'Experimental normgrcyt');
            hold off;
        
            title(['Comparison for Cell ID: ', num2str(cell_ids(i))]);
            xlabel('Time');
            ylabel('Concentration');
            legend show;
            grid on;
        end
end

%% ODE system for GR_alpha + GR_beta
% Define the ODE system
function dydt = crs_odes(t, y)

    % Define parameters
    kcn0 = 0.005;
    kcn1 = 0.02;
    MDex = 5.0;
    Dex0 = 100.0;
    gDex = 0.003;
    
    ba1 = 14e-5;
    da1 = 1e-5;
    ka1 = 0.01;
    da2 = 1e-6;
    
    bb1 = 14e-5;
    db1 = 1e-5;
    kb1 = 0.01;
    kb2 = 0.01;
    db2 = 1e-6;
    
    % Define input signal
    IDex = @(t) Dex0 * exp(-gDex * t);

    % Unpack variables
    x1 = y(1); % cytGR_a
    x2 = y(2); % nucGR_a
    x3 = y(3); % cytGR_b
    x4 = y(4); % nucGR_b

    % Reaction rates
    w1 = (kcn0 + kcn1 * IDex(t) / (MDex + IDex(t))) * x1;
    w2 = ba1;
    w3 = da1 * x1;
    w4 = ka1 * x2;
    w5 = da2 * x2;

    w6 = bb1;
    w7 = db1 * x3;
    w8 = kb1 * x3;
    w9 = kb2 * x4;
    w10 = db2 * x4;

    % ODEs
    dx1_dt = -w1 + w2 - w3 + w4;
    dx2_dt = w1 - w4 - w5;
    dx3_dt = w6 - w7 - w8 + w9;
    dx4_dt = w8 - w9 - w10;

    % Return derivatives
    dydt = [dx1_dt; dx2_dt; dx3_dt; dx4_dt];
end