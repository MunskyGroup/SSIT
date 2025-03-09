%% Fit Eric Ron's GR model to 3 concentrations of Dex: 1,10,100nM

close all 
clear
addpath(genpath('../../src'));
addpath('tmpPropensityFunctions');

    fitOptions = optimset('Display','iter','MaxIter',200);

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

    
    [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve;
    [FSPGrSoln,ModelGR.fspOptions.bounds] = ModelGR.solve(FSPGrSoln.stateSpace);

    % Define GR parameters
    GRpars = cell2mat(ModelGR.parameters(5:12,2))';  

    %% Associate GR Data with Different Instances of Model (10,100nm Dex)
    GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
        '10','10',102,'GR Fit (10nM Dex)';...
        '100','100',103,'GR Fit (100nM Dex)'};
    ModelGRparameterMap = cell(1,size(GRfitCases,1));
    ModelGRfit = cell(1,size(GRfitCases,1));

    % OLD GR DATA
    % ModelGRfit{1} = ModelGR.loadData("EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins_Dex_Conc_1.csv",...
    %            {'nucGR','normgrnuc';'cytGR','normgrcyt'});
    % ModelGRfit{2} = ModelGR.loadData("EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins_Dex_Conc_10.csv",...
    %            {'nucGR','normgrnuc';'cytGR','normgrcyt'});
    % ModelGRfit{3} = ModelGR.loadData("EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins_Dex_Conc_100.csv",...
    %            {'nucGR','normgrnuc';'cytGR','normgrcyt'});
    % disp(ModelGRfit{1}.dataSet.DATA(1:5,:))
    % disp(ModelGRfit{2}.dataSet.DATA(1:5,:))
    % disp(ModelGRfit{3}.dataSet.DATA(1:5,:))

    % NEW GR DATA
    % ModelGRfit{1} = ModelGR.loadData("EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03_dex_conc_1.csv",...
    %           {'nucGR','normGRnuc';'cytGR','normGRcyt'});
    % ModelGRfit{2} = ModelGR.loadData("EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03_dex_conc_10.csv",...
    %           {'nucGR','normGRnuc';'cytGR','normGRcyt'});
    % ModelGRfit{3} = ModelGR.loadData("EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03_dex_conc_100.csv",...
    %           {'nucGR','normGRnuc';'cytGR','normGRcyt'});
    % disp(ModelGRfit{1}.dataSet.DATA(1:5,:))
    % disp(ModelGRfit{2}.dataSet.DATA(1:5,:))
    % disp(ModelGRfit{3}.dataSet.DATA(1:5,:))
    
    for i=1:3
        ModelGRfit{i} = ModelGR.loadData("EricData/Gated_dataframe_Ron_020224_NormalizedGR_bins.csv",...
                                        {'nucGR','normgrnuc';'cytGR','normgrcyt'},...
                                        {'Dex_Conc',GRfitCases{i,2}});  

        disp(ModelGRfit{i}.dataSet.DATA(1:5,:))

        %ModelGRfit{i} = ModelGR.loadData("EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv",{'nucGR','normGRnuc';'cytGR','normGRcyt'},
        
        ModelGRfit{i}.parameters(13,:) = {'Dex0', str2num(GRfitCases{i,1})};
        ModelGRparameterMap(i) = {(1:8)};
    end

    disp(ModelGRfit{1}.dataSet.DATA(1:5,:))
    
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
         GRpars = combinedGRModel.maximizeLikelihood(GRpars, fitOptions);
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

%%
% dex_conc_values = [1, 10, 100];
% 
% for i = 1:3
%     % Load the full dataset first
%     fullData = readtable("EricData/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv");
% 
%     % Filter based on dex_conc value
%     filteredData = fullData(fullData.dex_conc == dex_conc_values(i), :);
% 
%     % Save filtered data to a temporary file (optional if loadData can take a table directly)
%     tempFile = sprintf("temp_filtered_data_%d.csv", dex_conc_values(i));
%     writetable(filteredData, tempFile);
% 
%     % Load the filtered data
%     combinedGRModel.SSITModels{i} = combinedGRModel.SSITModels{i}.loadData(tempFile, ...
%                 {'nucGR','normGRnuc';'cytGR','normGRcyt'});
% end


%%     STEP 1.F. -- Make Plots of GR Fit Results
makeGRPlots(combinedGRModel,GRpars)
