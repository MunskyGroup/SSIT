function fitModel2Data(app,fit_data_name,fit_results_file_name,Merged)
%% This function fits a model to data loaded in the data loading
% and fitting tab of the GUI. The user can fit the data using one of
% the three fitting methods listed in a drop down menu on the GUI.
arguments
    app
    fit_data_name = [];
    fit_results_file_name = [];
    Merged = false;
end
error('This function obsolete and set for deletion')

if Merged
    ssit.parest.fitMergedModel2Data(app,fit_data_name,fit_results_file_name);
    return
end

if nargin<=2
    if isempty(app.FspConstraintTable.Data)
        makeDefaultConstraints(app);
        readConstraintsForAdaptiveFsp(app);
    end
    
    %% Run initial FSP to get expansion.
    if isempty(app.pdo_parameters_table.Data)
        vars_to_fit1 = strcmp(app.fit_parameters_table.Data(:,3),'y');
        vars_to_fit2 = [];
        x0 = log10(cell2mat(app.fit_parameters_table.Data(vars_to_fit1,2)));
    else
        vars_to_fit1 = strcmp(app.fit_parameters_table.Data(:,3),'y');
        vars_to_fit2 = strcmp(app.pdo_parameters_table.Data(:,3),'y');
        x0 = [log10(cell2mat(app.fit_parameters_table.Data(vars_to_fit1,2)));...
            log10(cell2mat(app.pdo_parameters_table.Data(vars_to_fit2,2)))];
    end
    x0 = max(x0,-6);
    startingError = get_fit_error(10.^x0,app,vars_to_fit1,vars_to_fit2);
    disp(['Starting fit with initial log-likelihood of ',num2str(startingError)]);
    
    %% Copy necessary fields of app into structure for use in fitting.
    fields = {'ParEstFitTimesList.Items',...
        'ParEstFitTimesList.Value',...
        'FspInitCondField.Value',...
        'ModelParameterTable.Data',...app.ModelParameterTable.Data;
        'ModelInputTable.Data',...app.ModelInputTable.Data;
        'ReactionsTabOutputs.parameters',...app.ReactionsTabOutputs.parameters;
        'ReactionsTabOutputs.propensities',...app.ReactionsTabOutputs.propensities;
        'FspPrintTimesField.Value',...app.FspPrintTimesField.Value;
        'ReactionsTabOutputs.inputs',...app.ReactionsTabOutputs.inputs;
        'FspPiecewiseCheckBox',...app.FspPiecewiseCheckBox;
        'ReactionsTabOutputs.propensities',...app.ReactionsTabOutputs.propensities;        'FspUseMexCheckBox.Value',...app.FspUseMexCheckBox.Value;
        'FspConstraintTable.Data',...app.FspConstraintTable.Data;
        'DataLoadingAndFittingTabOutputs.boundIndex',...app.DataLoadingAndFittingTabOutputs.boundIndex;
        'FspTabOutputs.fConstraints',...app.FspTabOutputs.fConstraints;
        'FspTabOutputs.bounds',...app.FspTabOutputs.bounds;
        'ReactionsTabOutputs.stoichMatrix',...app.ReactionsTabOutputs.stoichMatrix;
        'DataLoadingAndFittingTabOutputs.dataTensor',...app.DataLoadingAndFittingTabOutputs.dataTensor;
        'DataLoadingAndFittingTabOutputs.fittingOptions',...app.DataLoadingAndFittingTabOutputs.fittingOptions;
        'fit_parameters_table.Data',...app.fit_parameters_table.Data;
        'FittingAlgorithmDropDown.Value',...app.FittingAlgorithmDropDown.Value;
        'FspIntegratorTolerance.Value',... app.FspIntegratorTolerance.Value;
        'DataLoadingAndFittingTabOutputs.fitOptions',...app.DataLoadingAndFittingTabOutputs.fitOptions;
        'DataLoadingAndFittingTabOutputs.priorOptions',...app.DataLoadingAndFittingTabOutputs.priorOptions;
        'DataLoadingAndFittingTabOutputs.J_LogLk',...app.DataLoadingAndFittingTabOutputs.J_LogLk;
        'PriorTypeDropDown.Value',...
        'FIMTabOutputs.distortionOperator',...
        'DistortionTypeDropDown.Value',...
        'FIMTabOutputs.PDOProperties',...
        'pdo_parameters_table.Data',...
        'ReactionsTabOutputs.varNames',...
        'FspTabOutputs.stateSpace',...
        'SpeciesForFitPlot',...
        'initApproxSS'};
    
    for i=1:length(fields)
        eval(['Loc_app_str.',fields{i},'=app.',fields{i},';']);
    end
    
    
    if min(eval(app.FspPrintTimesField.Value))<min(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)
        Loc_app_str.FspPrintTimesField.Value = ['[',num2str(min(eval(app.FspPrintTimesField.Value))),',',num2str(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times),']'];
    else
        Loc_app_str.FspPrintTimesField.Value = ['[',num2str(app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times),']'];
    end
    
    if app.SuppressFSPExpansionfasterbutmaybelessaccurateCheckBox.Value
        Loc_app_str.FspErrorTolField.Value = inf;
    else
        Loc_app_str.FspErrorTolField.Value = app.FspErrorTolField.Value;
    end
    
    Loc_app_str = struct(Loc_app_str);
    
    vars_to_fit = strcmp(Loc_app_str.fit_parameters_table.Data(:,3),'y');
    
    if nargin==2
        save(fit_data_name,'Loc_app_str','vars_to_fit','fields')
        return
    end
elseif nargin>=3
    fit_data_name = ['',fit_data_name,''];
    fit_results_file_name = ['',fit_results_file_name,''];
    load(fit_data_name,'Loc_app_str','vars_to_fit','fields');
end


if isempty(Loc_app_str.pdo_parameters_table.Data)
    vars_to_fit = strcmp(Loc_app_str.fit_parameters_table.Data(:,3),'y');
    x0 = log10(cell2mat(Loc_app_str.fit_parameters_table.Data(vars_to_fit,2)));
else
    vars_to_fit1 = strcmp(Loc_app_str.fit_parameters_table.Data(:,3),'y');
    vars_to_fit2 = strcmp(Loc_app_str.pdo_parameters_table.Data(:,3),'y');
    x0 = [log10(cell2mat(Loc_app_str.fit_parameters_table.Data(vars_to_fit1,2)));...
        log10(cell2mat(Loc_app_str.pdo_parameters_table.Data(vars_to_fit2,2)))];
end
x0 = max(x0,-6);
OBJLogLikelihood = @(x)-get_fit_error(10.^x,Loc_app_str,vars_to_fit1,vars_to_fit2);

switch Loc_app_str.PriorTypeDropDown.Value
    case 'None'
        OBJ = OBJLogLikelihood;
    case 'Normal'
        props = Loc_app_str.DataLoadingAndFittingTabOutputs.priorOptions.props;
        J = find(strcmp(Loc_app_str.fit_parameters_table.Data(:,3),'y'));
        for i=1:length(J)
            try
                mn(i) = eval(props.(['MEAN_',Loc_app_str.fit_parameters_table.Data{i,1}]));
                dev(i) = eval(props.(['STDV_',Loc_app_str.fit_parameters_table.Data{i,1}]));
            catch
                mn(i) = props.(['MEAN_',Loc_app_str.fit_parameters_table.Data{i,1}]);
                dev(i) = props.(['STDV_',Loc_app_str.fit_parameters_table.Data{i,1}]);
            end
        end
        if ~isempty(Loc_app_str.pdo_parameters_table.Data)
            J = find(strcmp(Loc_app_str.pdo_parameters_table.Data(:,3),'y'));
            for k=1:length(J)
                try
                    mn(i+k) = eval(props.(['MEAN_',Loc_app_str.pdo_parameters_table.Data{k,1}]));
                    dev(i+k) = eval(props.(['STDV_',Loc_app_str.pdo_parameters_table.Data{k,1}]));
                catch
                    mn(i+k) = props.(['MEAN_',Loc_app_str.pdo_parameters_table.Data{k,1}]);
                    dev(i+k) = props.(['STDV_',Loc_app_str.pdo_parameters_table.Data{k,1}]);
                end
            end
        end       
        OBJ = @(x)(OBJLogLikelihood(x) + sum((10.^x-mn(:)).^2./(2*dev(:).^2)));
    case 'LogNormal'
        props = Loc_app_str.DataLoadingAndFittingTabOutputs.priorOptions.props;
        J = find(strcmp(Loc_app_str.fit_parameters_table.Data(:,3),'y'));
        for i=1:length(J)
            try
                mn(i) = eval(props.(['LogMEAN_',Loc_app_str.fit_parameters_table.Data{i,1}]));
                dev(i) = eval(props.(['LogSTDV_',Loc_app_str.fit_parameters_table.Data{i,1}]));
            catch
                mn(i) = props.(['LogMEAN_',Loc_app_str.fit_parameters_table.Data{i,1}]);
                dev(i) = props.(['LogSTDV_',Loc_app_str.fit_parameters_table.Data{i,1}]);
            end
        end
        if ~isempty(Loc_app_str.pdo_parameters_table.Data)
            J = find(strcmp(Loc_app_str.pdo_parameters_table.Data(:,3),'y'));
            for k=1:length(J)
                try
                    mn(i+k) = eval(props.(['LogMEAN_',Loc_app_str.fit_parameters_table.Data{k,1}]));
                    dev(i+k) = eval(props.(['LogSTDV_',Loc_app_str.fit_parameters_table.Data{k,1}]));
                catch
                    mn(i+k) = props.(['LogMEAN_',Loc_app_str.fit_parameters_table.Data{k,1}]);
                    dev(i+k) = props.(['LogSTDV_',Loc_app_str.fit_parameters_table.Data{k,1}]);
                end
            end
        end
        OBJ = @(x)(OBJLogLikelihood(x) + sum((x-mn(:)).^2./(2*dev(:).^2)));
end

if isempty(Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions)
    Loc_app_str = ssit.parest.changeFitType(Loc_app_str);
end
switch Loc_app_str.FittingAlgorithmDropDown.Value
    case 'fminsearch'
        %         Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props
        xbest = fminsearch(OBJ,x0,Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props);
    case 'simulated annealing'
        rng('shuffle')
        xbest = ssit.parest.anneal(OBJ, x0, Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props);

    case 'genetic algorithm'
        rng('shuffle')
        OBJga = @(x)OBJ(x');
        LB = -5*ones(size(x0'));
        UB = 5*ones(size(x0'));
        Mut_Fun = @(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation)Mutation(parents,options,nvars,FitnessFcn,state,thisScore,thisPopulation,LB,UB);
        Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props.MutationFcn = Mut_Fun;
        X = repmat(x0',Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props.PopulationSize-1,1);
        Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props.InitialPopulationMatrix = [x0';X.*(1+0.1*randn(size(X)))];
        xbest = ga(OBJga,length(x0),[],[],[],[],LB,UB,[],...
            Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props);
        
    case 'particle swarm'
        rng('shuffle')
        OBJps = @(x)OBJ(x');
        LB = -5*ones(size(x0'));
        UB = 5*ones(size(x0'));
        initSwarm = repmat(x0',Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props.SwarmSize-1,1);
        initSwarm = [x0';initSwarm.*(1+0.1*randn(size(initSwarm)))];
        Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props.InitialSwarmMatrix = initSwarm;
                      
        xbest = particleswarm(OBJps,length(x0),LB,UB,...
            Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props);
        
    case 'Metropolis Hastings'
        rng('shuffle')
        OBJmh = @(x)-OBJ(x');
        options = Loc_app_str.DataLoadingAndFittingTabOutputs.fitOptions.props;
        if options.numChains==1

            if isstring(options.proposalDistribution)||ischar(options.proposalDistribution)
                options.proposalDistribution=eval(options.proposalDistribution);
            end

            [Loc_app_str.DataLoadingAndFittingTabOutputs.mhSamples,...
                Loc_app_str.DataLoadingAndFittingTabOutputs.mhAcceptance,...
                Loc_app_str.DataLoadingAndFittingTabOutputs.mhValue,xbest] = ...
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
            Loc_app_str.DataLoadingAndFittingTabOutputs.mhSamples = tmpMHSamp;
            Loc_app_str.DataLoadingAndFittingTabOutputs.mhAcceptance = tmpMHAcceptance;
            Loc_app_str.DataLoadingAndFittingTabOutputs.mhValue = tmpMHValue;
            clear tmpMH*
            
        end
        
        xbest = xbest';
        
    case 'Hamiltonian MC'
        Loc_app_str.FspErrorTolField.Value = inf;
        disp('Repressing FSP expansion during HMC');
        OBJmh = @(x)-OBJ(x);
        % Create hamiltinian MC sampler for model and data.
        hmc = hmcSampler(OBJmh,x0,'UseNumericalGradient',true);%,'StepSize',0.01);
        
        % Estimate the MLE
        % [x0] = estimateMAP(hmc,'VerbosityLevel',1,);
        
        % Tune a HMC sampler
        warning('off','MATLAB:nearlySingularMatrix')
        [smp] = tuneSampler(hmc,'Start',x0,'VerbosityLevel',2);
                
        [Loc_app_str.DataLoadingAndFittingTabOutputs.hMCSamples,...
            Loc_app_str.DataLoadingAndFittingTabOutputs.hMCEndpoint,...
            Loc_app_str.DataLoadingAndFittingTabOutputs.hMCAccRatio] = ...
            drawSamples(smp,'VerbosityLevel',2,...
            'StartPoint',x0);
        
end

[~,Loc_app_str] = get_fit_error(10.^xbest,Loc_app_str,vars_to_fit1,vars_to_fit2);

%% Put results of the fits back into the app.
if nargin==1
    
    app.FspTabOutputs.solutions = Loc_app_str.FspTabOutputs.solutions;
    app.FspTabOutputs.solutions = Loc_app_str.FspTabOutputs.solutions;
    app.DataLoadingAndFittingTabOutputs.fitResults = Loc_app_str.DataLoadingAndFittingTabOutputs.fitResults;
    app.fit_parameters_table.Data = Loc_app_str.fit_parameters_table.Data;
    app.ModelParameterTable.Data(:,2) = app.fit_parameters_table.Data(:,2);
    
    app.FspPrintTimesField.Value = Loc_app_str.FspPrintTimesField.Value ;
    
    try
        app.DataLoadingAndFittingTabOutputs.mhSamples = Loc_app_str.DataLoadingAndFittingTabOutputs.mhSamples;
        app.DataLoadingAndFittingTabOutputs.mhValue = Loc_app_str.DataLoadingAndFittingTabOutputs.mhValue;
    catch
    end
    
    %     runFsp(app);
    %     makePlotOfData(app)
    
    ssit.parest.updateModelSolveAndCompareToData(app);
    makeSeparatePlotOfData(app);
    figure(app.UIFigure);

elseif nargin>=3
    save(fit_results_file_name,'Loc_app_str','vars_to_fit','fields')
end
end

function [fit_error,app] = get_fit_error(x,app,vars_to_fit1,vars_to_fit2)
arguments
    x
    app
    vars_to_fit1
    vars_to_fit2=[];
end
app.fit_parameters_table.Data(vars_to_fit1,2) = num2cell(x(1:sum(vars_to_fit1)));
if ~isempty(vars_to_fit2)
    app.pdo_parameters_table.Data(vars_to_fit2,2) = num2cell(x(sum(vars_to_fit1)+1:end));
end
app.ReactionsTabOutputs.parameters(:,2) = app.fit_parameters_table.Data(:,2);
app = ssit.parest.updateModelSolveAndCompareToData(app);
fit_error = app.DataLoadingAndFittingTabOutputs.J_LogLk;
end

function [fit_log_likelihood,fit_log_likelihood_gradient] = get_fit_error_and_sensitivity(x,app,vars_to_fit1,vars_to_fit2)
arguments
    x
    app
    vars_to_fit1
    vars_to_fit2=[];
end
app.fit_parameters_table.Data(vars_to_fit1,2) = num2cell(x(1:sum(vars_to_fit1)));
if ~isempty(vars_to_fit2)
    app.pdo_parameters_table.Data(vars_to_fit2,2) = num2cell(x(sum(vars_to_fit1)+1:end));
end
app.ReactionsTabOutputs.parameters(:,2) = app.fit_parameters_table.Data(:,2);
app = ssit.parest.updateModelSolveAndCompareToData(app,true);
fit_log_likelihood = app.DataLoadingAndFittingTabOutputs.J_LogLk;
fit_log_likelihood_gradient = app.DataLoadingAndFittingTabOutputs.dLogLk_dpar;
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
