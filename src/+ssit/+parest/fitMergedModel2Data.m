function fitMergedModel2Data(app,SSITapp,fit_data_name,fit_results_file_name,isBackground)
%% This function fits a model to data loaded in the data loading
% and fitting tab of the GUI. The user can fit the data using one of
% the three fitting methods listed in a drop down menu on the GUI.
arguments
    app
    SSITapp
    fit_data_name = [];
    fit_results_file_name = [];
    isBackground = false;
end
if isempty(fit_data_name)||(isBackground&&isempty(fit_results_file_name))
    %% Copy necessary fields of app into structure for use in fitting.
    MergeModSettings.fit_parameters_table.Data = app.fit_parameters_table.Data;
    MergeModSettings.DataLoadingAndFittingTabOutputs.ModelMerge = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge;    
    MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions = app.DataLoadingAndFittingTabOutputs.fitOptions;
    MergeModSettings.FittingAlgorithmDropDown.Value = app.FittingAlgorithmDropDown.Value;
    MergeModSettings = struct(MergeModSettings);
    
    vars_to_fit = strcmp(MergeModSettings.fit_parameters_table.Data(:,3),'y');
    MergeModSettings.OBJ = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.OverallObjective;

    if isBackground&&isempty(fit_results_file_name)
        save([fit_data_name],'MergeModSettings','vars_to_fit')
        return
    end
elseif ~isempty(fit_data_name)  
    load(fit_data_name,'MergeModSettings','vars_to_fit');
end

OBJ = MergeModSettings.OBJ;
x0 = log10(cell2mat(MergeModSettings.fit_parameters_table.Data(vars_to_fit,2)));
x0 = max(x0,-6);

% Run initial parameter values to determine FSP bounds.
for i=1:length(SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FitResults)
    FitResults{i} = SSITapp.DataLoadingAndFittingTabOutputs.ModelMerge.FitResults{i}(x0);
    fspBounds{i} = FitResults{i}{4};
end
% Reset Merge model with current bounds.
ssit.parest.ModelMergeMerge(app,SSITapp,false,fspBounds,false)

if isempty(MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions)
    MergeModSettings = ssit.parest.changeFitType(MergeModSettings);
end
switch MergeModSettings.FittingAlgorithmDropDown.Value
    case 'fminsearch'
        xbest = fminsearch(OBJ,x0,MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props);
    case 'simulated annealing'
        rng('shuffle')
        xbest = ssit.parest.anneal(OBJ, x0, MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props);
    case 'genetic algorithm'
        rng('shuffle')
        OBJga = @(x)OBJ(x');
        LB = -5*ones(size(x0'));
        UB = 5*ones(size(x0'));
        Mut_Fun = @(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation)Mutation(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,LB,UB);
        MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props.MutationFcn = Mut_Fun;
        X = repmat(x0',MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props.PopulationSize-1,1);
        MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props.InitialPopulationMatrix = [x0';X.*(1+0.1*randn(size(X)))];
        xbest = ga(OBJga,length(x0),[],[],[],[],LB,UB,[],...
            MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props);

    case 'particle swarm'
        rng('shuffle')
        OBJps = @(x)OBJ(x');
        LB = -5*ones(size(x0'));
        UB = 5*ones(size(x0'));
        initSwarm = repmat(x0',MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props.SwarmSize-1,1);
        initSwarm = [x0';initSwarm.*(1+0.1*randn(size(initSwarm)))];
        MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props.InitialSwarmMatrix = initSwarm;
        optionsA = MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props;
       
        fldNames = fieldnames(optionsA);
        options = optimoptions('particleSwarm');
        for i=1:length(fldNames)
            if ~isempty(optionsA.(fldNames{i}))
                options.(fldNames{i}) = optionsA.(fldNames{i});
            end
        end
        xbest = particleswarm(OBJps,length(x0),LB,UB,options);
              
    case 'Metropolis Hastings'
        rng('shuffle')
        OBJmh = @(x)-OBJ(x');
        options = MergeModSettings.DataLoadingAndFittingTabOutputs.fitOptions.props;
        if options.numChains==1
            [MergeModSettings.DataLoadingAndFittingTabOutputs.mhSamples,...
                MergeModSettings.DataLoadingAndFittingTabOutputs.mhAcceptance,...
                MergeModSettings.DataLoadingAndFittingTabOutputs.mhValue,xbest] = ...
                ssit.parest.metropolisHastingsSample(x0',options.numberOfSamples,...
                'logpdf',OBJmh,'proprnd',options.proposalDistribution,'symmetric',options.isPropDistSymmetric,...
                'thin',options.thin,'nchain',1,'burnin',options.burnIn,...
                'progress',options.progress);
        else
            try
                parpool
            catch
            end
            options.progress=0;
            clear tmpMH*
            parfor iChain = 1:options.numChains
                [mhSamples, mhAcceptance, mhValue,xbest,fbest] = ...
                    ssit.parest.metropolisHastingsSample(x0',options.numberOfSamples,...
                    'logpdf',OBJmh,'proprnd',options.proposalDistribution,'symmetric',options.isPropDistSymmetric,...
                    'thin',options.thin,'nchain',1,'burnin',options.burnIn,...
                    'progress',options.progress);
                tmpMHSamp(iChain) = {mhSamples};
                tmpMHAcceptance(iChain) = {mhAcceptance};
                tmpMHValue(iChain) = {mhValue};
                tmpMHxbest(iChain) = {xbest};
                tmpMHfbest(iChain) = fbest;
            end
            [~,jBest] = max(tmpMHfbest);
            xbest = tmpMHxbest{jBest};
            MergeModSettings.DataLoadingAndFittingTabOutputs.mhSamples = tmpMHSamp;
            MergeModSettings.DataLoadingAndFittingTabOutputs.mhAcceptance = tmpMHAcceptance;
            MergeModSettings.DataLoadingAndFittingTabOutputs.mhValue = tmpMHValue;
            clear tmpMH*
            
        end
        
        xbest = xbest';
        
    case 'Hamiltonian MC'
        MergeModSettings.FspErrorTolField.Value = inf;
        disp('Repressing FSP expansion during HMC');
        OBJmh = @(x)-OBJ(x);
        hmc = hmcSampler(OBJmh,x0,'UseNumericalGradient',true,'StepSize',0.01);
        [MAPpars] = estimateMAP(hmc,'VerbosityLevel',1);
        [smp] = tuneSampler(hmc,'Start',MAPpars);
        [MergeModSettings.DataLoadingAndFittingTabOutputs.hMCSamples,...
            MergeModSettings.DataLoadingAndFittingTabOutputs.hMCEndpoint,...
            MergeModSettings.DataLoadingAndFittingTabOutputs.hMCAccRatio] = ...
            drawSamples(smp,'VerbosityLevel',2,...
            'StartPoint',x0);      
end

% MergeModSettings = get_fit_error(10.^xbest,MergeModSettings,vars_to_fit);
MergeModSettings.fit_parameters_table.Data(vars_to_fit,2) = num2cell(10.^xbest);

%% Put results of the fits back into the app.
if isempty(fit_data_name)
    app.fit_parameters_table.Data = MergeModSettings.fit_parameters_table.Data;
    try
        app.DataLoadingAndFittingTabOutputs.mhSamples = MergeModSettings.DataLoadingAndFittingTabOutputs.mhSamples;
        app.DataLoadingAndFittingTabOutputs.mhValue = MergeModSettings.DataLoadingAndFittingTabOutputs.mhValue;
    catch
    end
elseif ~isempty(fit_results_file_name)
    save(fit_results_file_name,'MergeModSettings')
end

if isBackground
    ssit.parest.SaveProject(MergeModSettings,MergeModSettings,[],fit_results_file_name);
else
    ssit.parest.SaveProject(app,SSITapp,[],fit_results_file_name);
end

end

function mut_Chil = Mutation(parents,~,~,~,~,~,this_Pop,LB,UB)
%% Mutation function for use in the genetic algorithm fits.
mut_Chil = this_Pop(parents,:).*...
    (1+(randn(size(this_Pop(parents,:)))>0.5).*...
    randn(size(this_Pop(parents,:)))/10^(1+3*rand))+...
    10^(-6*rand)*randn*(rand(size(this_Pop(parents,:)))>0.8);

for i=1:size(mut_Chil,1)
    mut_Chil(i,:) = max([mut_Chil(i,:);LB]);
    mut_Chil(i,:) = min([mut_Chil(i,:);UB]);
end
end
