classdef SSIT

    properties
        parameters = {'k',10; 'g',0.2};   % List of parameters and their values.
        species = {'x1'}; % List of species to be used in model (x1,x2,...)
        stoichiometry = [1,-1]; % Stoichiometry matrix
        propensityFunctions = {'k'; 'g*x1'} % List of proensity functions
        inputExpressions = {}; % List of time varying input signals (I1,I2,...)
        customConstraintFuns = {}; % User suppled constraint functions for FSP.
        fspOptions = struct('fspTol',0.001,'fspIntegratorRelTol',1e-2,...
            'fspIntegratorAbsTol',1e-4, 'odeSolver','auto', 'verbose',false,...
            'bounds',[],'usePiecewiseFSP',false,...
            'initApproxSS',false,...
            'escapeSinks',[],...
            'constantJacobian',false, ...
            'constantJacobianTime',1.1); % Options for FSP solver.
        sensOptions = struct('solutionMethod','forward','useParallel',true);
        % Options for FSP-Sensitivity solver.
        ssaOptions = struct('Nexp',1,'nSimsPerExpt',100,'useTimeVar',false,...
            'signalUpdateRate',[],'useParalel',false,...
            'verbose',false); % Options for SSA solver
        pdoOptions = struct('unobservedSpecies',[],'PDO',[]);
        % Options for FIM analyses
        fittingOptions = struct('modelVarsToFit','all','pdoVarsToFit',[],...
            'timesToFit','all','logPrior',[])
        initialCondition = [0]; % Initial condition for species [x1;x2;...]
        initialProbs = 1; % Probability mass of states given in init. cond.
        initialTime = 0;
        tSpan = linspace(0,10,21); % Times at which to find solutions
        solutionScheme = 'FSP' % Chosen solutuon scheme ('FSP','SSA')
        modelReductionOptions = struct('useModReduction',false,'reductionType','None') % Settings for
        % model reduction tools.
        dataSet = [];
        useHybrid = false
        hybridOptions = struct('upstreamODEs',[]);
        propensitiesGeneral = [];% Processed propensity functions for use in solvers

    end

    properties (Dependent)
        fspConstraints % FSP Constraint Functions
        pars_container  % Container for parameters
    end

    methods
        function obj = SSIT(modelFile)
            % SSIT - create an instance of the SSIT class.
            % Arguments:
            %   modelFile (optional) -- create  from specified template:
            %               {'Empty',
            %                'BirthDeath',    % 1 species example
            %                'CentralDogma',  % 2 species example
            %                'ToggleSwitch',  % 2 species example
            %                'Repressilator', % 3 species example
            %                'BurstingSpatialCentralDogma'}    % 4 species example
            % Example:
            %   F = SSIT('CentralDogma'); % Generate model for
            %                           %transcription and translation.
            arguments
                modelFile = [];
            end
            % SSIT Construct an instance of the SSIT class
            addpath(genpath('../src'));
            if ~isempty(modelFile)
                obj = pregenModel(obj,modelFile);
            end
        end

        function Pars_container = get.pars_container(obj)
            if ~isempty(obj.parameters)
                Pars_container = containers.Map(obj.parameters(:,1), obj.parameters(:,2));
            else
                Pars_container =[];
            end
        end

       
        function obj = formPropensitiesGeneral(obj,prefixName)
            arguments
                obj
                prefixName = [];
            end
            % This function starts the process to write m-file for each
            % propensity function.
            % if ~strcmp(obj.solutionScheme,'SSA')
                if strcmp(obj.solutionScheme,'ode')
                    obj.useHybrid = true;
                    obj.hybridOptions.upstreamODEs = obj.species;
                end
                n_reactions = length(obj.propensityFunctions);
                % Propensity for hybrid models will include
                % solutions from the upstream ODEs.
                sm = cell(1,n_reactions);
                logicTerms = cell(1,n_reactions);
                logCounter = 0;
                for i = 1:n_reactions
                    st = obj.propensityFunctions{i};
                    for jI = 1:size(obj.inputExpressions,1)
                        st = strrep(st,obj.inputExpressions{jI,1},['(',obj.inputExpressions{jI,2},')']);
                    end
                    [st,logicTerms{i},logCounter] = ssit.Propensity.stripLogicals(st,obj.species,logCounter);
                    sm{i} = str2sym(st);
                end

                if obj.useHybrid
                    PropensitiesGeneral = ssit.Propensity.createAsHybridVec(sm, obj.stoichiometry,...
                        obj.parameters, obj.species, obj.hybridOptions.upstreamODEs, logicTerms, prefixName);
                else
                    PropensitiesGeneral = ssit.Propensity.createAsHybridVec(sm, obj.stoichiometry,...
                        obj.parameters, obj.species, [], logicTerms, prefixName);
                end

                obj.propensitiesGeneral = PropensitiesGeneral;

            % elseif strcmp(obj.solutionScheme,'SSA')
            %     PropensitiesGeneral = ssit.SrnModel.processPropensityStrings(obj.propensityFunctions,...
            %         obj.inputExpressions,...
            %         obj.pars_container,...
            %         propenType,...
            %         obj.species);
            % end

        end
%%
        function constraints = get.fspConstraints(obj)
            % Makes a list of FSP constraints that can be used by the FSP
            % solver.
            if obj.useHybrid
                stochasticSpecies = setdiff(obj.species,obj.hybridOptions.upstreamODEs);
            else
                stochasticSpecies = obj.species;
            end

            nSpecies = length(stochasticSpecies);
            Data = cell(nSpecies*2,3);
            for i = 1:nSpecies
                Data(i,:) = {['-x',num2str(i)],'<',0};
                Data(nSpecies+i,:) = {['x',num2str(i)],'<',1};
            end
            for i = 1:length(obj.customConstraintFuns)
                Data(2*nSpecies+i,:) = {obj.customConstraintFuns{i},'<',1};
            end
            constraints.f = readConstraintsForAdaptiveFsp([], stochasticSpecies, Data);
            if isempty(obj.fspOptions.bounds)||size(Data,1)~=length(obj.fspOptions.bounds)
                constraints.b = [Data{:,3}]';
                obj.fspOptions.bounds = constraints.b;
            else
                constraints.b = obj.fspOptions.bounds;
            end

            if ~isempty(obj.fspOptions.escapeSinks)
                nEscape = length(obj.fspOptions.escapeSinks.f);
                escapeData = cell(nEscape,3);
                for i = 1:nEscape
                    escapeData(i,:) = {obj.fspOptions.escapeSinks.f{i},'<',1};
                end
                constraints.fEscape = readConstraintsForAdaptiveFsp([], stochasticSpecies, escapeData);
                constraints.bEscape = obj.fspOptions.escapeSinks.b;
            else
                constraints.fEscape = [];
                constraints.bEscape = [];
            end
        end

        %% Model Building Functions
        function [obj] = pregenModel(obj,modelFile)
            % pregenModel - creates a pregenerated model from a template:
            % Possible Templates include:
            %   Empty -- nothing
            %   BirthDeath -- one species 'x1' with birth rate 'k' and
            %       death rate 'g'
            %   CentralDogma -- Time varying 2-species model with:
            %       mRNA species 'x1' with birth rate 'kr*I(t)' and
            %       degradation rate 'gr'. Protein species 'x2' with
            %       translation rate 'kr' and degradation rate 'gp'.
            %   ToggleSwitch -- two proteins that prepress one another with
            %       non-linear functions.
            switch modelFile
                case 'Empty'
                    obj.parameters = {};
                    obj.species = {};
                    obj.stoichiometry = [];
                    obj.propensityFunctions = {};
                    obj.initialCondition = [];
                case 'BirthDeath'
                    obj.parameters = {'k',10;'g',0.2};
                    obj.species = {'x1'};
                    obj.stoichiometry = [1,-1];
                    obj.propensityFunctions = {'k';'g*x1'};
                    obj.initialCondition = [0];
                case 'CentralDogma'
                    obj.parameters = {'kr',10;'gr',1;'kp',1;'gp',0.1};
                    obj.species = {'x1';'x2'};
                    obj.stoichiometry = [1,-1,0, 0;...
                        0, 0,1,-1];
                    obj.propensityFunctions = {'kr';'gr*x1';'kp*x1';'gp*x2'};
                    obj.initialCondition = [0;0];
                case 'BurstingGene'
                    obj.parameters = {'kon',1;'koff',1;'kr',1;'gr',0.1};
                    obj.species = {'x1';'x2'};
                    obj.stoichiometry = [1,-1,0, 0;...
                        0, 0,1,-1];
                    obj.propensityFunctions = {'kon*(1-x1)';'koff*x1';'kr*x1';'gr*x2'};
                    obj.initialCondition = [0;0];
                case 'CentralDogmaTV'
                    obj.parameters = {'kr',10;'gr',1;'kp',1;'gp',0.1;'omega',2*pi/5};
                    obj.species = {'x1';'x2'};
                    obj.stoichiometry = [1,-1,0, 0;...
                        0, 0,1,-1];
                    obj.propensityFunctions = {'kr';'gr*I1*x1';'kp*x1';'gp*x2'};
                    obj.initialCondition = [0;0];
                    obj.inputExpressions = {'I1','1+cos(omega*t)'};
                case 'ToggleSwitch'
                    obj.parameters = {'kb',10;'ka',80;'M',20;'g',1};
                    obj.species = {'x1';'x2'};
                    obj.stoichiometry = [1,-1,0, 0;...
                        0, 0,1,-1];
                    obj.propensityFunctions = {'kb+ka*M^3/(M^3+x2^3)';...
                        'g*x1';...
                        'kb+ka*M^3/(M^3+x1^3)';...
                        'g*x2'};
                    obj.initialCondition = [0;0];
                    obj.customConstraintFuns = {'(x1-3).^2.*(x2-3).^2'};
                case 'ToggleSwitch2'
                    obj.parameters = {'ka1',4;'kb1',80;'kd1',1;'k1',20;...
                        'ka2',4;'kb2',80;'kd2',1;'k2',20;...
                        'ket',0.1;'ks',1;'kg',1};
                    obj.species = {'x1';'x2'};
                    obj.stoichiometry = [1,-1,0, 0;...
                        0, 0,1,-1];
                    obj.propensityFunctions = {'ket*(ka1+((kb1*(k1^3))/((k1^3)+(x2)^3)))';...
                        '(kd1+((ks*kg)/(1+ks)))*(x1)';...
                        'ket*(ka2+((kb2*(k2^3))/((k2^3)+(x1)^3)))';...
                        'kd2*(x2)'};
                    obj.initialCondition = [0;0];
                    obj.customConstraintFuns = {'(x1-3).^2.*(x2-3).^2'};
                case 'Repressilator'
                    obj.parameters = {'kn0',0;'kn1',25;'a',5;'n',6;'g',1};
                    obj.species = {'x1';'x2';'x3'};
                    obj.stoichiometry = [1,0,0,-1,0,0;...
                        0,1,0,0,-1,0;...
                        0,0,1,0,0,-1];
                    obj.propensityFunctions = {'kn0+kn1*(1/(1+a*(x2^n)))';...
                        'kn0+kn1*(1/(1+a*(x3^n)))';...
                        'kn0+kn1*(1/(1+a*(x1^n)))';...
                        'g*x1';...
                        'g*x2';...
                        'g*x3'};
                    obj.initialCondition = [30;0;0];
                    obj.customConstraintFuns = {'(x1-3).^2.*(x2-3).^2*(x3-3).^2'};

                case 'RepressilatorGenes'
                    obj.parameters = {'kn0',0;'kn1',25;'kb',2000;'ku',10;'g',1};
                    obj.species = {'x1';'x2';'x3';'x4';'x5';'x6';'x7';'x8';'x9'};
                    obj.stoichiometry = zeros(9,12);
                    obj.stoichiometry(1,1:2) = [-1 1];
                    obj.stoichiometry(2,1:2) = [1 -1];
                    obj.stoichiometry(6,1:2) = [-3 3];
                    obj.stoichiometry(3,3) =  1;
                    obj.stoichiometry(3,4) = -1;
                    obj.propensityFunctions(1:4) = {'kb*x1*x6*(x6-1)/2*(x6-2)/6';'ku*x2';'kn0*x2+kn1*x1';'g*x3'};
                    obj.stoichiometry(4,5:6) = [-1 1];
                    obj.stoichiometry(5,5:6) = [1 -1];
                    obj.stoichiometry(9,5:6) = [-3 3];
                    obj.stoichiometry(6,7) =  1;
                    obj.stoichiometry(6,8) = -1;
                    obj.propensityFunctions(5:8) = {'kb*x4*x9*(x9-1)/2*(x9-2)/6';'ku*x5';'kn0*x5+kn1*x4';'g*x6'};
                    obj.stoichiometry(7,9:10) = [-1 1];
                    obj.stoichiometry(8,9:10) = [1 -1];
                    obj.stoichiometry(3,9:10) = [-3 3];
                    obj.stoichiometry(9,11) =  1;
                    obj.stoichiometry(9,12) = -1;
                    obj.propensityFunctions(9:12) = {'kb*x7*x3*(x3-1)/2*(x3-2)/6';'ku*x8';'kn0*x8+kn1*x7';'g*x9'};
                    obj.initialCondition = [1;0;30;0;1;0;0;1;0];
                    obj.customConstraintFuns = {'(x3-3).^3.*(x6-3).^3.*(x9-3).^3'};

                case 'BurstingSpatialCentralDogma'
                    obj.parameters = {'kon',1;'koff',2;...
                        'kr',5;'grn',0.1;'kt',0.5;...
                        'grc',0.1;...
                        'kp',1;'gp',0.1};
                    obj.species = {'x1';'x2';'x3';'x4'};
                    obj.stoichiometry = [1,-1,0,0,0,0,0,0;...
                        0,0,1,-1,-1,0,0,0;...
                        0,0,0,0,1,-1,0,0;...
                        0,0,0,0,0,0,1,-1];
                    obj.propensityFunctions = {'kon*(1-x1)';'koff*x1';...
                        'kr*x1';'grn*x2';'kt*x2';...
                        'grc*x3';...
                        'kp*x3';'gp*x4'};
                    obj.initialCondition = [0;0;0;0];
                    obj.customConstraintFuns = {};

            end
            obj.propensitiesGeneral = [];
        end

        function [obj] = createModelFromSBML(obj,sbmlFile,scaleVolume)
             arguments
                obj
                sbmlFile
                scaleVolume = false
            end
           % This function allows one to create a model directly from an
            % SBML file.
            % Example:
            %      Model = SSIT();
            %      Model = Model.createModelFromSBML('../SBML_test_cases/00010/00010-sbml-l1v2.xml');
            %      [fspSoln] = Model.solve;
            %      Model.makePlot(fspSoln,'meansAndDevs')
            sbmlobj = sbmlimport(sbmlFile);
            [obj] = createModelFromSimBiol(obj,sbmlobj,scaleVolume);
            obj.propensitiesGeneral = [];

        end

        function [obj] = createModelFromSimBiol(obj,sbmlobj,scaleVolume)
           % This function allows one to create a model directly from an
            % simBiology object.
            % Example:
            %      Model = SSIT();
            %      Model = Model.createModelFromSimBiol(sbmlobj);
            %      [fspSoln] = Model.solve;
            %      Model.makePlot(fspSoln,'meansAndDevs')
            arguments
                obj
                sbmlobj
                scaleVolume = false
            end
            nR = length(sbmlobj.Reactions);
            nS = length(sbmlobj.Species);

            % Extract species names and stoichiometry
            [obj.stoichiometry, obj.species] = getstoichmatrix(sbmlobj);

            % Extract parameter names
            nP = length(sbmlobj.Parameter);
            obj.parameters = {};
            for i = 1:nP
                obj.parameters{i,1} = sbmlobj.Parameter(i).Name;
                obj.parameters{i,2} = sbmlobj.Parameter(i).Value;
            end

            if length(sbmlobj.Compartments)>1
                error('SSIT Tools not yet set up to support multi-compartment models.')
            end

            obj.propensityFunctions={};
            for i = 1:nR
                reactRate = sbmlobj.Reactions(i).ReactionRate;
                reactRate = strrep(reactRate,'compartment*','');
                if contains(reactRate,'time')
                    obj.propensityFunctions{i,1} = strrep(reactRate,'time','Ig');
                    obj.inputExpressions = {'Ig','t'};
                else
                    obj.propensityFunctions{i,1} = reactRate;
                end

            end

            if scaleVolume
                % Replace species numbers (Xi) with concentrations (Xi/Volume).
                for i = 1:nR
                    for j = 1:nS
                        obj.propensityFunctions{i,1} = strrep(obj.propensityFunctions{i,1},...
                            obj.species{j},['(',obj.species{j},'/Volume)']);
                    end
                    obj.propensityFunctions{i,1} = [obj.propensityFunctions{i,1},'*Volume'];
                end

                % Scale Initial Condition and Volume to remove fractional
                % concentrations.
                frac = false;
                scl = 0;
                for i = 1:nS
                    if rem(sbmlobj.Species(i).Value,1)~=0
                        frac = true;
                    end
                    scl = max(scl,sbmlobj.Species(i).Value);
                end
                if frac
                    scl = round(100/scl);
                    disp(['Fractional species values detected.  Scaling by Vol=',num2str(scl),' and rounding.'])
                    disp(' ')
                    IC(1:nS,1) = round(scl*[sbmlobj.Species.Value]);
                end
                obj.parameters(end+1,:) = {'Volume',scl};
            else
                IC(1:nS,1) = [sbmlobj.Species.Value];
            end
            obj.initialCondition = IC;
            obj.summarizeModel;
            obj.propensitiesGeneral = [];

        end

        function exportToSBML(obj,modelName)
            % This function exports the model to an SBML file called
            % <modelName>.
            arguments
                obj
                modelName
            end
            sbModel = exportSimBiol(obj);
            sbmlexport(sbModel, modelName)
        end

        function sbModel = exportSimBiol(obj,verifyAndPlot)
            % This function converts the model to a simple simbiology model
            % and returns that simbiology object. 
            % Arguments:
            %   verifyAndPlot (true/false0) -- option to verify the model
            %       and run simBiology to make a plot of its results.
            %
            % Outputs:
            %   smModel -- the resulting simBiology model.
            arguments
                obj
                verifyAndPlot = false;
            end

            sbModel = sbiomodel('simpleModel');

            s = cell(1,length(obj.species));        
            for is = 1:size(obj.stoichiometry,1)
                s{is} = addspecies(sbModel,obj.species{is},obj.initialCondition(is));
            end

            for is = 1:size(obj.parameters,1)
                p{is} = addparameter(sbModel,obj.parameters{is,1},obj.parameters{is,2});
            end

            % Parse time varying components in the reaction rate equations. 
            props = obj.propensityFunctions;
            for is = 1:size(obj.inputExpressions,1)
                tvComp = strrep(obj.inputExpressions{is,2},'t','time');
                for ir = 1:length(props)
                    props{ir} = strrep(props{ir},obj.inputExpressions{is,1},['(',tvComp,')']);
                end
            end
           
            for ir=1:size(obj.stoichiometry,2)
                strReactants = [];
                strProducts =  [];
                for is = 1:size(obj.stoichiometry,1)
                    if obj.stoichiometry(is,ir)<0
                        strReactants =[strReactants,'+ ',[num2str(-obj.stoichiometry(is,ir))],' ',obj.species{is}];
                    elseif obj.stoichiometry(is,ir)>0
                        strProducts =[strProducts,'+ ',[num2str(obj.stoichiometry(is,ir))],' ',obj.species{is}];
                    end
                end
                if isempty(strProducts); strProducts = '  null '; end
                if isempty(strReactants); strReactants = '  null '; end
                rxn = [strReactants(3:end),' -> ',strProducts(3:end)];
                RXN{is} = addreaction(sbModel,rxn,'ReactionRate',props{ir});
            end

            if verifyAndPlot
                verify(sbModel)
                csObj = getconfigset(sbModel,'active');
                set(csObj,'Stoptime',max(obj.tSpan));
                [t,x,names] = sbiosimulate(sbModel);
                plot(t,x);
                xlabel('Time');
                ylabel('Amount');
                legend(names);
            end
        end

        function [obj] = addSpecies(obj,newSpecies,initialCond)
            % addSpecies - add new species to reaction model.
            % example:
            %     F = SSIT;
            %     F = F.addSpecies('x2');
            arguments
                obj
                newSpecies
                initialCond = [];
            end
            obj.species =  [obj.species;newSpecies];
            obj.stoichiometry(end+1,:) = 0;
            if isempty(initialCond)
                initialCond = zeros(size(newSpecies,1),1);
            end
            obj.initialCondition = [obj.initialCondition;initialCond];
            obj.propensitiesGeneral = [];

        end

        function [obj] = addParameter(obj,newParameters)
            % addParameter - add new parameter to reaction model
            % example:
            %     F = SSIT;
            %     F = F.addParameter({'kr',0.1})
            obj.parameters =  [obj.parameters;newParameters];
        end

        function [obj] = addReaction(obj,newPropensity,newStoichVector)
            % addParameter - add new reaction to reaction model
            % example:
            %     F = SSIT;
            %     F = F.addReaction({'kr*x1'},[0;-1]);  % Add reaction
            %     x1->x1+x2 with rate kr.
            obj.propensityFunctions =  [obj.propensityFunctions;newPropensity];
            obj.stoichiometry =  [obj.stoichiometry,newStoichVector];
            obj.propensitiesGeneral = [];
        end

        function [obj] = calibratePDO(obj,dataFileName,measuredSpecies,...
                trueColumns,measuredColumns,pdoType,showPlot,parGuess)
            % This function calibrates the PDO to match some provided data
            % for 'true' and 'observed' spot numbers.
            arguments
                obj
                dataFileName
                measuredSpecies
                trueColumns
                measuredColumns
                pdoType = 'AffinePoiss'
                showPlot = false
                parGuess=[];
            end

            obj.pdoOptions.type = pdoType;
            %             app.DistortionTypeDropDown.Value = obj.pdoOptions.type;
            %             app.FIMTabOutputs.PDOProperties.props = obj.pdoOptions.props;

            Tab = readtable(dataFileName);
            dataNames = Tab.Properties.VariableNames;
            DATA = table2cell(Tab);

            if isempty(parGuess)
                lambdaTemplate = obj.findPdoError(pdoType);
            else
                lambdaTemplate=parGuess;
            end

            lambda = [];
            maxSize = zeros(1,length(obj.species));
            options = optimset('display','none');
            for i=1:length(obj.species)
                if sum(strcmp(measuredSpecies,obj.species{i}))==1
                    k = find(strcmp(measuredSpecies,obj.species{i}));
                    jTrue = find(strcmp(dataNames,trueColumns{k}));
                    jObsv = find(strcmp(dataNames,measuredColumns{k}));
                    xTrue = [DATA{:,jTrue}]';
                    xObsv = [DATA{:,jObsv}]';
                    maxSize(i)=max(xTrue);
                    objPDO = @(x)-obj.findPdoError(pdoType,x,xTrue,xObsv);
                    lambdaNew = fminsearch(objPDO,lambdaTemplate,options);
                    if showPlot
                        [~,PDO] = obj.findPdoError(pdoType,lambdaNew,xTrue,xObsv);
                        Z = max(-200,log10(PDO));
                        figure
                        contourf([0:size(PDO,2)-1],[0:size(PDO,1)-1],Z);
                        colorbar
                        hold on
                        scatter(xTrue,xObsv,100,'sk','filled')
                        set(gca,'fontsize',15)
                        legend('PDO','Data')
                    end
                else
                    maxSize(i)=0;
                    lambdaNew = 0*lambdaTemplate;
                end
                lambda = [lambda,lambdaNew];
            end
            obj.pdoOptions.props.PDOpars = lambda;
            obj.pdoOptions.PDO = obj.generatePDO(obj.pdoOptions,lambda,[],[],maxSize);
        end
        
        function [pdo] = generatePDO(obj,pdoOptions,paramsPDO,fspSoln,variablePDO,maxSize)
            arguments
                obj
                pdoOptions
                paramsPDO = []
                fspSoln = []
                variablePDO =[]
                maxSize=[];
            end
            app.DistortionTypeDropDown.Value = pdoOptions.type;
            app.FIMTabOutputs.PDOProperties.props = pdoOptions.props;
            % Separate into observed and unobserved species.
            Nd = length(obj.species);
            indsUnobserved=[];
            indsObserved=[];
            for i=1:Nd
                if ~isempty(obj.pdoOptions.unobservedSpecies)&&max(contains(obj.pdoOptions.unobservedSpecies,obj.species{i}))
                    indsUnobserved=[indsUnobserved,i];
                else
                    indsObserved=[indsObserved,i];
                end
            end
            [~,pdo] = ssit.pdo.generatePDO(app,paramsPDO,fspSoln,indsObserved,variablePDO,maxSize);
        end

        function [logL,P] = findPdoError(obj,pdoType,lambda,True,Distorted)
            % This function calculates the likelihood of observed data
            % given true data and an assumed PDO model.
            arguments
                obj
                pdoType = 'AffinePoiss';
                lambda = [];
                True = [];
                Distorted =[];
            end

            if nargin<=2
                switch pdoType
                    case 'AffinePoiss'
                        logL = [1 5 0.5];
                end
                return
            end

            % Computes likelihood of observed data given the model of affine poisson
            % extra spot counting and probability of measurmeent failure.
            NmaxTrue = max(True);
            NmaxObs = max(Distorted);

            switch pdoType
                case 'AffinePoiss'
                    Np = ceil(max(NmaxObs,lambda(2)+lambda(3)*NmaxTrue));
            end
            P = zeros(Np+1,NmaxTrue+1);

            for xi = 0:NmaxTrue
                switch pdoType
                    case 'AffinePoiss'
                        P(1:Np+1,xi+1) = pdf('poiss',[0:Np]',max(lambda(1),lambda(2)+lambda(3)*xi));
                end
            end

            % compute likelihood of observed given true
            logP = max(log(P),-100);
            logL = 0;
            for i = 1:length(True)
                logL = logL + logP(Distorted(i)+1,True(i)+1);
            end

            % apply constraints
            switch pdoType
                case 'AffinePoiss'
                    logL = logL-1e4*(lambda(1)<0);
            end
        end

        function summarizeModel(obj)
            arguments
                obj;
            end
            % Show the model species
            nS = size(obj.stoichiometry,1);
            disp('Species:')
            for i = 1:nS
                if ~isempty(obj.hybridOptions.upstreamODEs)&&max(contains(obj.hybridOptions.upstreamODEs,obj.species{i}))
                    disp(['     ',obj.species{i},'; IC = ',num2str(obj.initialCondition(i)),';  upstream ODE']);
                else
                    disp(['     ',obj.species{i},'; IC = ',num2str(obj.initialCondition(i)),';  discrete stochastic']);
                end
            end
            disp(' ')

            % Show the model stoichiometries and propensity functions
            disp('Reactions:')
            nR = size(obj.stoichiometry,2);
            for iR = 1:nR
                s = obj.stoichiometry(:,iR);
                disp(['  Reaction ',num2str(iR),':'])
                jReactant = find(s<0);
                jProd = find(s>0);
                if isempty(jReactant)
                    reactTxt = 'NULL';
                else
                    reactTxt = [num2str(-s(jReactant(1))),'*',obj.species{jReactant(1)}];
                    for i = 2:length(jReactant)
                        reactTxt = [reactTxt,' + ',num2str(-s(jReactant(i))),'*',obj.species{jReactant(i)}];
                    end
                end
                if isempty(jProd)
                    prodTxt = 'NULL';
                else
                    prodTxt = [num2str(s(jProd(1))),'*',obj.species{jProd(1)}];
                    for i = 2:length(jProd)
                        prodTxt = [prodTxt,' + ',num2str(s(jProd(i))),'*',obj.species{jProd(i)}];
                    end
                end
                disp(['     s',num2str(iR),': ',reactTxt, ' --> ', prodTxt])

                disp(['     w',num2str(iR),': ',obj.propensityFunctions{iR}])

            end

            if ~isempty(obj.inputExpressions)
                disp(' ')
                disp('Input Signals:')
                nI = size(obj.inputExpressions,1);
                for i = 1:nI
                    disp(['     ',obj.inputExpressions{i,1},'(t) = ',obj.inputExpressions{i,2}])
                end
            end

            disp(' ')
            disp('Model Parameters:')
            disp(obj.parameters)

        end

        %% Model Analysis Functions
        function [Solution, bConstraints] = solve(obj,stateSpace,saveFile,fspSoln)
            arguments
                obj
                stateSpace = [];
                saveFile=[];
                fspSoln=[];
            end
            % solve - solve the model using the specified method in
            %    obj.solutionScheme
            % Example:
            %   F = SSIT('ToggleSwitch')
            %   F.solutionScheme = 'FSP'
            %   [soln,bounds] = F.solve;  % Returns the solution and the
            %                             % bounds for the FSP projection
            %   F.solutionScheme = 'fspSens'
            %   [soln,bounds] = F.solve;  % Returns the sensitivity and the
            %                             % bounds for the FSP projection
            % See also: SSIT.makePlot for information on how to visualize
            % the solution data.
            if obj.initialTime>obj.tSpan(1)
                error('First time in tspan cannot be earlier than the initial time.')
            elseif obj.initialTime~=obj.tSpan(1)
                %                 warning('First time in tspan is not the same as initial time.')
                obj.tSpan = unique([obj.initialTime,obj.tSpan]);
            end

            if isempty(obj.propensitiesGeneral)
                obj = formPropensitiesGeneral(obj);
            end

            if obj.modelReductionOptions.useModReduction
                if ~isfield(obj.modelReductionOptions,'phi')
                    error('Model Reduction Matrices have not yet been Defined.')
                end
                useReducedModel = true;
                modRedTransformMatrices.phi = obj.modelReductionOptions.phi;
                modRedTransformMatrices.phi_inv = obj.modelReductionOptions.phi_inv;
                modRedTransformMatrices.phiScale = obj.modelReductionOptions.phiScale;
                modRedTransformMatrices.phiPlot = obj.modelReductionOptions.phiPlot;
            else
                useReducedModel = false;
                modRedTransformMatrices = [];
            end

            switch obj.solutionScheme
                case 'FSP'
                    if ~isempty(stateSpace)&&size(stateSpace.states,2)~=stateSpace.state2indMap.Count
                        error('HERE')
                    end

                    % specificPropensities = SSIT.parameterizePropensities(obj.propensitiesGeneral,[obj.parameters{:,2}]');

                    [Solution.fsp, bConstraints,Solution.stateSpace] = ssit.fsp.adaptiveFspSolve(obj.tSpan,...
                        obj.initialCondition,...
                        obj.initialProbs,...
                        obj.stoichiometry, ...
                        obj.propensitiesGeneral, ...
                        [obj.parameters{:,2}]', ...
                        obj.fspOptions.fspTol, ...
                        obj.fspConstraints.f, ...
                        obj.fspConstraints.b,...
                        obj.fspOptions.verbose, ...
                        obj.fspOptions.fspIntegratorRelTol, ...
                        obj.fspOptions.fspIntegratorAbsTol, ...
                        obj.fspOptions.odeSolver,stateSpace,...
                        obj.fspOptions.usePiecewiseFSP,...
                        obj.fspOptions.initApproxSS,...
                        obj.species,...
                        useReducedModel,modRedTransformMatrices, ...
                        obj.useHybrid,obj.hybridOptions,...
                        obj.fspConstraints.fEscape,obj.fspConstraints.bEscape, ...
                        obj.fspOptions.constantJacobian);

                case 'SSA'
                    Solution.T_array = obj.tSpan;
                    Nt = length(Solution.T_array);
                    nSims = obj.ssaOptions.Nexp*obj.ssaOptions.nSimsPerExpt*Nt;
                    W = obj.propensitiesGeneral;
                    if obj.ssaOptions.useParalel
                        trajs = zeros(length(obj.species),...
                            length(obj.tSpan),nSims);% Creates an empty Trajectories matrix from the size of the time array and number of simulations
                        parfor isim = 1:nSims
                            if obj.ssaOptions.verbose
                                disp(['completed sim number: ',num2str(isim)])
                            end
                            trajs(:,:,isim) = ssit.ssa.runSingleSsa(obj.initialCondition,...
                                obj.stoichiometry,...
                                W,...
                                obj.tSpan,...
                                obj.ssaOptions.useTimeVar,...
                                obj.ssaOptions.signalUpdateRate,...
                                [obj.parameters{:,2}]');
                        end
                        Solution.trajs = trajs;
                    else
                        Solution.trajs = zeros(length(obj.species),...
                            length(obj.tSpan),nSims);% Creates an empty Trajectories matrix from the size of the time array and number of simulations
                        for isim = 1:nSims
                            Solution.trajs(:,:,isim) = ssit.ssa.runSingleSsa(obj.initialCondition,...
                                obj.stoichiometry,...
                                W,...
                                obj.tSpan,...
                                obj.ssaOptions.useTimeVar,...
                                obj.ssaOptions.signalUpdateRate,...
                                [obj.parameters{:,2}]');
                        end
                    end
                    disp([num2str(nSims),' SSA Runs Completed'])
                    try
                        if ~isempty(obj.pdoOptions.PDO)
                            Solution.trajsDistorted = zeros(length(obj.species),...
                                length(obj.tSpan),nSims);% Creates an empty Trajectories matrix from the size of the time array and number of simulations
                            for iS = 1:length(obj.species)
                                PDO = obj.pdoOptions.PDO.conditionalPmfs{iS};
                                nDpossible = size(PDO,1);
                                Q = Solution.trajs(iS,:,:);
                                for iD = 1:length(Q(:))
                                    Q(iD) = randsample([0:nDpossible-1],1,true,PDO(:,Q(iD)+1));
                                end
                                Solution.trajsDistorted(iS,:,:) = Q;
                            end
                            disp('PDO applied to SSA results')
                        end
                        if ~isempty(saveFile)
                            A = table;
                            for j=1:Nt
                                A.time((j-1)*obj.ssaOptions.nSimsPerExpt+1:j*obj.ssaOptions.nSimsPerExpt) = obj.tSpan(j);
                                for i = 1:obj.ssaOptions.Nexp
                                    for k=1:obj.ssaOptions.nSimsPerExpt
                                        for s = 1:size(Solution.trajs,1)
                                            warning('off')
                                            A.(['exp',num2str(i),'_s',num2str(s)])((j-1)*obj.ssaOptions.nSimsPerExpt+k) = ...
                                                Solution.trajs(s,j,(i-1)*Nt*obj.ssaOptions.nSimsPerExpt+(j-1)*obj.ssaOptions.nSimsPerExpt+k);
                                            if ~isempty(obj.pdoOptions.PDO)
                                                A.(['exp',num2str(i),'_s',num2str(s),'_Distorted'])((j-1)*obj.ssaOptions.nSimsPerExpt+k) = ...
                                                    Solution.trajsDistorted(s,j,(i-1)*Nt*obj.ssaOptions.nSimsPerExpt+(j-1)*obj.ssaOptions.nSimsPerExpt+k);
                                            end
                                        end
                                    end
                                end
                            end
                            writetable(A,saveFile)
                            disp(['SSA Results saved to ',saveFile])
                        end
                    catch
                        pause;
                    end
                case 'fspSens'
                    if ~isempty(obj.parameters)
                        model = ssit.SrnModel(obj.stoichiometry,...
                            obj.propensityFunctions,...
                            obj.parameters(:,1),...
                            obj.inputExpressions);
                        app.ReactionsTabOutputs.parameters = obj.parameters(:,1);
                    else
                        model = ssit.SrnModel(obj.stoichiometry,...
                            obj.propensityFunctions,...
                            [],...
                            obj.inputExpressions);
                        app.ReactionsTabOutputs.parameters = [];
                    end
                    app.ReactionsTabOutputs.varNames = obj.species;
                    [Solution.sens, bConstraints] = ...
                        ssit.sensitivity.computeSensitivity(model,...
                        obj.parameters,...
                        obj.propensitiesGeneral,...
                        obj.tSpan,...
                        obj.fspOptions.fspTol,...
                        obj.initialCondition,...
                        1.0,...
                        obj.fspConstraints.f,...
                        obj.fspConstraints.b,...
                        [], obj.fspOptions.verbose, 0,...
                        obj.sensOptions.solutionMethod,...
                        app,stateSpace,...
                        obj.fspOptions.usePiecewiseFSP,...
                        obj.fspOptions.initApproxSS,...
                        obj.species,...
                        obj.sensOptions.useParallel,...
                        fspSoln,...
                        useReducedModel,modRedTransformMatrices, ...
                        obj.useHybrid,obj.hybridOptions,...
                        obj.fspConstraints.fEscape,obj.fspConstraints.bEscape);
                    %                     app.SensFspTabOutputs.solutions = Solution.sens;
                    %                     app.SensPrintTimesEditField.Value = mat2str(obj.tSpan);
                    %                     Solution.plotable = exportSensResults(app);

                case 'ode'
                    
                    % specificPropensities = SSIT.parameterizePropensities(obj.propensitiesGeneral,[obj.parameters{:,2}]');
                    
                    [~,Solution.ode] = ssit.moments.solveOde2(obj.initialCondition, obj.tSpan, ...
                        obj.stoichiometry, obj.propensitiesGeneral,  [obj.parameters{:,2}]', obj.fspOptions.initApproxSS);
            end
        end

        function sampleDataFromFSP(obj,fspSoln,saveFile)
            Solution.T_array = obj.tSpan;
            Nt = length(Solution.T_array);
            nSims = obj.ssaOptions.nSimsPerExpt*obj.ssaOptions.Nexp;
            Solution.trajs = zeros(length(obj.species),...
                length(obj.tSpan),nSims);% Creates an empty Trajectories matrix
            % from the size of the time array and number of simulations
            for it = 1:length(obj.tSpan)
                clear PP
                PP = double(fspSoln.fsp{it}.p.data);
                clear w
                w(:) = PP(:); w(w<0)=0;
                [I1,I2,I3,I4,I5] =  ind2sub(size(PP),randsample(length(w), nSims, true, w ));
                for iSp = 1:length(obj.species)
                    eval(['Solution.trajs(iSp,it,:) = I',num2str(iSp),'-1;']);
                end
            end
            if ~isempty(obj.pdoOptions.PDO)
                Solution.trajsDistorted = zeros(length(obj.species),...
                    length(obj.tSpan),nSims);% Creates an empty Trajectories matrix from the size of the time array and number of simulations
                for iS = 1:length(obj.species)
                    PDO = obj.pdoOptions.PDO.conditionalPmfs{iS};
                    nDpossible = size(PDO,1);
                    Q = Solution.trajs(iS,:,:);
                    for iD = 1:length(Q(:))
                        Q(iD) = randsample([0:nDpossible-1],1,true,PDO(:,Q(iD)+1));
                    end
                    Solution.trajsDistorted(iS,:,:) = Q;
                end
                disp('PDO applied to FSP Samples')
            end
            if ~isempty(saveFile)
                A = table;
                for it=1:Nt
                    A.time((it-1)*obj.ssaOptions.nSimsPerExpt+1:it*obj.ssaOptions.nSimsPerExpt) = obj.tSpan(it);
                    for ie = 1:obj.ssaOptions.Nexp
                        %                          for is=1:obj.ssaOptions.nSimsPerExpt
                        for s = 1:size(Solution.trajs,1)
                            warning('off')
                            A.(['exp',num2str(ie),'_s',num2str(s)])((it-1)*obj.ssaOptions.nSimsPerExpt+(1:obj.ssaOptions.nSimsPerExpt)) = ...
                                Solution.trajs(s,it,(ie-1)*obj.ssaOptions.nSimsPerExpt+(1:obj.ssaOptions.nSimsPerExpt));
                            if ~isempty(obj.pdoOptions.PDO)
                                A.(['exp',num2str(ie),'_s',num2str(s),'_Distorted'])((it-1)*obj.ssaOptions.nSimsPerExpt+(1:obj.ssaOptions.nSimsPerExpt)) = ...
                                    Solution.trajsDistorted(s,it,(ie-1)*obj.ssaOptions.nSimsPerExpt+(1:obj.ssaOptions.nSimsPerExpt));
                            end
                        end
                        %                          end
                    end
                end
                writetable(A,saveFile)
                disp(['FSP Samples saved to ',saveFile])
            end
        end

        function [fimResults,sensSoln] = computeFIM(obj,sensSoln,scale,MHSamples)
            % computeFIM - computes FIM at all time points.
            % Arguments:
            %   sensSoln - (optional) previously compute FSP Sensitivity.
            %              Automatically computed if not provided.
            %   scale - ('lin'(default) or 'log') Choice of FIM based on
            %            linear parameters or their natural logarithm
            %   MHSamples - (optional) set of parameter sets at which to calculate
            %           the FIM.
            % Outputs:
            %   fimResults - FIM at each time point in obj.tSpan
            %   sensSoln - FSP Sensitivity.
            arguments
                obj
                sensSoln = [];
                scale = 'lin';
                MHSamples = [];
            end

            if ~isempty(MHSamples)
                % For FIM calculation 
                nSamples = size(MHSamples,1);
                Nt = length(obj.tSpan);
                fimResults = cell(Nt,nSamples);
                if isempty(sensSoln)||length(sensSoln)~=nSamples
                    sensSoln = cell(1,nSamples);
                end
                if strcmp(obj.fittingOptions.modelVarsToFit,'all')
                    obj.fittingOptions.modelVarsToFit = 1:size(obj.parameters,1);
                end
                if nargout == 2
                    saveSens = true;
                else
                    saveSens = false;
                end

                for i=1:nSamples
                    objTMP = obj;
                    objTMP.parameters(objTMP.fittingOptions.modelVarsToFit,2) = ...
                        num2cell(MHSamples(i,:));
                    if saveSens
                        [fimResults(:,i),sensSoln{i}] = computeFIM(obj,sensSoln{i},scale);
                    else
                        fimResults(:,i) = objTMP.computeFIM(sensSoln{i},scale);
                    end
                end
            else

                if isempty(sensSoln)
                    disp({'Running Sensitivity Calculation';'You can skip this step by providing sensSoln.'})
                    obj.solutionScheme = 'fspSens';
                    [sensSoln] = obj.solve;
                    sensSoln = sensSoln.sens;
                end

                % Separate into observed and unobserved species.
                Nd = length(obj.species);
                indsUnobserved=[];
                indsObserved=[];
                for i=1:Nd
                    if ~isempty(obj.pdoOptions.unobservedSpecies)&&max(contains(obj.pdoOptions.unobservedSpecies,obj.species{i}))
                        indsUnobserved=[indsUnobserved,i];
                    else
                        indsObserved=[indsObserved,i];
                    end
                end

                % compute FIM for each time point
                fimResults = {};
                for it=length(sensSoln.data):-1:1
                    if isempty(indsUnobserved)
                        F = ssit.fim.computeSingleCellFim(sensSoln.data{it}.p, sensSoln.data{it}.S, obj.pdoOptions.PDO);
                    else
                        % Remove unobservable species.
                        redS = sensSoln.data{it}.S;
                        for ir = 1:length(redS)
                            redS(ir) = sensSoln.data{it}.S(ir).sumOver(indsUnobserved);
                        end
                        F = ssit.fim.computeSingleCellFim(sensSoln.data{it}.p.sumOver(indsUnobserved), redS, obj.pdoOptions.PDO);
                    end
                    fimResults{it,1} = F;
                end

                if strcmp(scale,'log')
                    for it=length(sensSoln.data):-1:1
                        fimResults{it,1} = diag([obj.parameters{:,2}])*...
                            fimResults{it,1}*...
                            diag([obj.parameters{:,2}]);
                    end
                end

            end
        end

        function [fimTotal,mleCovEstimate,fimMetrics] = evaluateExperiment(obj,...
                fimResults,cellCounts)
            % This function evaluates the provided experiment design (in
            % "cellCounts" and produces an array of FIMs (one for each
            % parameter set.
            arguments
                obj
                fimResults
                cellCounts
            end
            Ns = size(fimResults,2);
            Nt = size(fimResults,1);
            Np = size(fimResults{1,1},1);
            fimTotal = cell(1,Ns);
            mleCovEstimate = cell(1,Ns);

            for is=1:Ns
                fimTotal{is} = 0*fimResults{1,is};

                for it=1:Nt
                    fimTotal{is} = fimTotal{is} + cellCounts(it)*fimResults{it,is};
                end

                if nargout>=2
                    % Estimate MLE covariance
                    if rank(fimTotal{is})<Np
                        disp(['FIM has rank ',num2str(rank(fimTotal{is})),' and is not invertable for this experiment design'])
                        mleCovEstimate{1,is} = NaN*ones(Np);
                    else
                        mleCovEstimate{1,is} = fimTotal{is}^-1;
                    end
                end
            end

            if nargout>=3
                for is = Ns:-1:1
                    % Compute FIM metrics.
                    fimMetrics.det(1,is) = det(fimTotal{is});
                    fimMetrics.trace(1,is) = trace(fimTotal{is});
                    fimMetrics.minEigVal(1,is) = min(eig(fimTotal{is}));
                end
            end
        end

        function [NcDNewDesign] = optimizeCellCounts(obj,fims,nCellsTotalNew,FIMMetric,NcGuess,NcFixed,NcMax)
            % This function optimizes the number of cells per time point
            % according to the user-provide metric. 
            % 
            % INPUTS:
            %
            % 'fims' is either [Nt x 1] cell array containing the FIM matrices
            % for each of the Nt time points, or it is a [Nt x Ns] cell
            % array containing the FIM for each combination of Nt time
            % points and Ns different parameter sets.
            %
            % 'nCellsTotalNew' is the total number of cells the user wishes to
            % measure, spead out among the Nt time points.
            % 
            % 'FIMmetric' is the type of optimization that the user
            % desires. Allowable metrics are:
            %   'Determinant' - maximize the expected determinant of the FIM
            %   'DetCovariance' - minimize the expected determinant of MLE covariance. 
            %   'Smallest Eigenvalue' - maximize the smallest e.val of the
            %       FIM
            %   'Trace' - maximize the trace of the FIM
            %   '[<i1>,<i2>,...]' - minimize the determinant of the inverse
            %       FIM for the specified indices. All other parameters are
            %       assumed to be free.
            %   'TR[<i1>,<i2>,...]' - maximize the determinant of the  FIM
            %       for the specified indices. Only the parameters in
            %       obj.fittingOptions.modelVarsToFit are assumed to be
            %       free.
            %
            % 'Nc' is an optimal guess for the optimal experiment design.
            %
            % 'NcFixed' is a minimal number of cells to measure at each
            %      time point.  This is useful for subsequent experiment
            %      design where you already have measured cells in the
            %
            % 'NcMax' maximum total number of cells allowed for each time
            %       point.  This is useful in simulated experiment design,
            %       where there are only so many cells available in the
            %       real data.
            %
            % OUTPUTS:
            % 'Nc' is the optimized experiment design (number of cells to
            % measure at each point in time.
            %
            arguments
                obj
                fims
                nCellsTotalNew
                FIMMetric = 'Smallest Eigenvalue';
                NcGuess = [];
                NcFixed = [];
                NcMax = []
            end
            switch FIMMetric
                case 'Determinant'
                    met = @(A)-det(A);
                case 'DetCovariance'
                    met = @(A)det(inv(A));
                case 'Smallest Eigenvalue'
                    met = @(A)-min(eig(A));
                case 'Trace'
                    met = @(A)-trace(A);
                otherwise
                    if strcmp(FIMMetric(1:2),'TR')
                        k = eval(FIMMetric(3:end));
                        met = @(A)det(inv(A(k,k)));
                    else  % all parameters are free.
                        k = eval(FIMMetric);
                        ek = zeros(length(k),length(fims{1}));
                        ek(1:length(k),k) = eye(length(k));
                        met = @(A)det(ek*inv(A)*ek');
                    end
            end
            NT = size(fims,1);
            NS = size(fims,2);

            if isempty(NcFixed)
                NcFixed = zeros(1,NT);
            end

            if isempty(NcMax)
                NcMax = inf*ones(1,NT);
            end
            
            if isempty(NcGuess)
                NcGuess = NcFixed;
                NcGuess(1)=NcGuess(1)+nCellsTotalNew;
            else
                NcGuess = NcFixed+NcGuess;
            end

            Converged = 0;
            while Converged==0
                Converged = 1;
                for i = 1:NT
                    while NcGuess(i)>NcFixed(i)
                        Ncp = NcGuess;
                        Ncp(i) = Ncp(i)-1;
                        k = SSIT.findBestMove(fims,Ncp,met,NcMax);
                        if k==i
                            break
                        end
                        NcGuess = Ncp;
                        NcGuess(k)=NcGuess(k)+1;
                        Converged = 0;
                    end
                end
            end
            NcDNewDesign = NcGuess - NcFixed;
        end

        %% Data Loading and Fitting
        function [obj] = loadData(obj,dataFileName,linkedSpecies,conditions)
            arguments
                obj
                dataFileName
                linkedSpecies
                conditions = {};
                % Data conditions that can be used to filter out data that
                % do not meet specifications.
                % Example:
                %     conditions = {'Rep_num','1'}  : Only the data in the
                %     'Rep_num' column that is exactly equal to '1' will be
                %     kept in the data set.
            end
            obj.dataSet =[];
            Tab = readtable(dataFileName);
            obj.dataSet.dataNames = Tab.Properties.VariableNames;
            obj.dataSet.DATA = table2cell(Tab);

            obj.dataSet.linkedSpecies = linkedSpecies;

            Q = contains(obj.dataSet.dataNames,{'time','Time','TIME'});
            if sum(Q)==1
                obj.dataSet.app.ParEstFitTimesList.Items = {};
                obj.dataSet.app.ParEstFitTimesList.Value = {};
                col_time = find(Q);
                obj.dataSet.app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_time_index = col_time;
                obj.dataSet.app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times = sort(unique(cell2mat(obj.dataSet.DATA(:,col_time))));
                for i=1:length(obj.dataSet.app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times)
                    obj.dataSet.app.ParEstFitTimesList.Items{i} = num2str(obj.dataSet.app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times(i));
                    obj.dataSet.app.ParEstFitTimesList.Value{i} = num2str(obj.dataSet.app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times(i));
                end
                % We need to make sure that the fitting times are included in the solution times.
            else
                error('Provided data set does not have required column named "time"')
            end
            obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTable = obj.dataSet.DATA;

            Nd = length(obj.species);
            nCol = length(obj.dataSet.dataNames);

            obj.dataSet.app.DataLoadingAndFittingTabOutputs.marginalMatrix = ...
                zeros(Nd+3,nCol);

            % auto-detect and record 'time' column
            Itime = Nd+1;
            Jtime = find(Q);
            obj.dataSet.app.DataLoadingAndFittingTabOutputs.marginalMatrix(Itime,Jtime) = 1;
%             obj.dataSet.times = unique([obj.dataSet.DATA{:,Jtime}]);

            % record linked species
            for i=1:size(linkedSpecies,1)
                J = find(strcmp(obj.dataSet.dataNames,linkedSpecies{i,2}));
                I = find(strcmp(obj.species,linkedSpecies{i,1}));
                obj.dataSet.app.DataLoadingAndFittingTabOutputs.marginalMatrix(I,J)=1;
            end

            % set up conditionals
            obj.dataSet.app.DataLoadingAndFittingTabOutputs.conditionOnArray = string(zeros(1,length(nCol)));
            for i=1:size(conditions,1)
                J = find(strcmp(obj.dataSet.dataNames,conditions{i,1}));
                obj.dataSet.app.DataLoadingAndFittingTabOutputs.conditionOnArray(J) = conditions{i,2};
            end

            % set to marginalize over everything else
            obj.dataSet.app.DataLoadingAndFittingTabOutputs.marginalMatrix(Nd+3,:) = ...
                sum(obj.dataSet.app.DataLoadingAndFittingTabOutputs.marginalMatrix)==0;

            obj.dataSet.app.SpeciesForFitPlot.Items = obj.species;
            [obj.dataSet.app,obj.dataSet.times] = filterAndMarginalize([],[],obj.dataSet.app);

%             obj.dataSet.times = unique([obj.dataSet.DATA{:,Jtime}]);

            sz = size(obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor);
            obj.dataSet.nCells=zeros(sz(1),1);
            for i=1:sz(1)
                if length(sz)==2
                    obj.dataSet.nCells(i) = sum(double(obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(i,:)),'all');
                elseif length(sz)==3
                    obj.dataSet.nCells(i) = sum(double(obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(i,:,:)),'all');
                elseif length(sz)==4
                    obj.dataSet.nCells(i) = sum(double(obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(i,:,:,:)),'all');
                elseif length(sz)==5
                    obj.dataSet.nCells(i) = sum(double(obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(i,:,:,:,:)),'all');
                end
            end

            obj.tSpan = unique([obj.initialTime,obj.dataSet.times]);

            % Calculate the means
            obj.dataSet.mean = zeros(sz(1),length(sz)-1);
            for i=1:sz(1)
                for j=2:length(sz)
                    tmpInt{j-1} = [1:sz(j)];
                end
                TMP = squeeze(double(obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(i,tmpInt{:})));
                for j = 1:length(sz)-1
                    if length(sz)>2
                        H = sum(TMP,[1:j-1,j+1:length(sz)-1]);
                    else
                        H = TMP;
                    end
                    obj.dataSet.mean(i,j) = sum([0:length(H)-1]'.*H(:))/sum(H);
                end
            end

        end

        function [logL,gradient] = minusLogL(obj,pars,stateSpace,computeSensitivity)
            [logL,gradient] = computeLikelihood(obj,exp(pars),stateSpace,computeSensitivity);
            logL = -logL;
            gradient = -gradient.*exp(pars);
        end

        function [logLode] = computeLikelihoodODE(obj,pars,SIG)
            arguments
                obj
                pars = [];
                SIG = [];
            end
        
            if strcmp(obj.fittingOptions.modelVarsToFit,'all')
                indsParsToFit = [1:length(obj.parameters)];
            else
                indsParsToFit = obj.fittingOptions.modelVarsToFit;
            end
            nModelPars = length(indsParsToFit);

            if isempty(pars)
                pars = [obj.parameters{indsParsToFit,2}];
            end

            if ~isempty(obj.fittingOptions.logPrior)
                logPrior = sum(obj.fittingOptions.logPrior(pars));
            else
                logPrior = 0;
            end

            originalPars = obj.parameters;
            obj.tSpan = unique([obj.initialTime,obj.tSpan,obj.dataSet.times]);
            [~,IA,~] = intersect(obj.tSpan,obj.dataSet.times);

            % Update Model and PDO parameters using supplied guess
            obj.parameters(indsParsToFit,2) =  num2cell(pars(1:nModelPars));

            obj.solutionScheme = 'ode'; % Chosen solutuon scheme
            for i=1:size(obj.parameters,1)
                obj.parameters{i,2} = round(obj.parameters{i,2},12);
            end
            solutions = obj.solve;  % Solve the ODE analysis

            obj.parameters =  originalPars;

            % Need to add likelihood calculation here.
            nt = length(IA);
%             ns = length(obj.species);

            for i = 1:size(obj.dataSet.linkedSpecies,1)
                J(i) = find(strcmp(obj.species,obj.dataSet.linkedSpecies{i,1}));
            end
            nds = length(J);

            if isempty(SIG)
                SIG = eye(nt*nds);
            end
            nc = repmat(obj.dataSet.nCells,nds,1);

            vm = zeros(nt*nds,1); 
            tmp = solutions.ode(IA,J);
            vm(:) = tmp(:);
            
            vd = zeros(nt*nds,1); vd(:) = obj.dataSet.mean(:);
            
            vm = real(vm);

            logLode = -1/2*(sqrt(nc)'.*(vd-vm)')*SIG^(-1)*((vd-vm).*sqrt(nc));
            logLode = logLode+logPrior;
        end

        function [logL,gradient,fitSolutions] = computeLikelihood(obj,pars,stateSpace,computeSensitivity)
            arguments
                obj
                pars = [];
                stateSpace =[];
                computeSensitivity = false;
            end

            if ~isempty(stateSpace)&&size(stateSpace.states,2)~=stateSpace.state2indMap.Count
                stateSpace =[];
            end

            if strcmp(obj.fittingOptions.modelVarsToFit,'all')
                indsParsToFit = [1:length(obj.parameters)];
            else
                indsParsToFit = obj.fittingOptions.modelVarsToFit;
            end
            nModelPars = length(indsParsToFit);

            if strcmp(obj.fittingOptions.pdoVarsToFit,'all')
                indsPdoParsToFit = [1:length(obj.pdoOptions.props.ParameterGuess)];
            else
                indsPdoParsToFit = obj.fittingOptions.pdoVarsToFit;
            end
            nPdoPars = length(indsPdoParsToFit);

            if isempty(pars)
                pars = [obj.parameters{indsParsToFit,2}];
            end

            if ~isempty(obj.fittingOptions.logPrior)
                logPrior = sum(obj.fittingOptions.logPrior(pars));
            else
                logPrior = 0;
            end

            originalPars = obj.parameters;
            %             obj.tSpan = unique([obj.initialTime,obj.dataSet.times]);
            obj.tSpan = unique([obj.initialTime,obj.tSpan]);

            % Update Model and PDO parameters using supplied guess
            obj.parameters(indsParsToFit,2) =  num2cell(pars(1:nModelPars));

            if computeSensitivity&&nargout>=2
                obj.solutionScheme = 'fspSens'; % Chosen solutuon scheme ('FSP','SSA')
                [solutions] = obj.solve(stateSpace);  % Solve the FSP analysis
            else
                obj.solutionScheme = 'FSP'; % Chosen solutuon scheme ('FSP','SSA')
                [solutions] = obj.solve(stateSpace);  % Solve the FSP analysis
            end
            obj.parameters =  originalPars;

            if ~isempty(pars)&&nPdoPars>0
                obj.pdoOptions.props.ParameterGuess(indsPdoParsToFit) = pars(nModelPars+1:end);
                obj.pdoOptions.PDO = obj.generatePDO(obj.pdoOptions,[],solutions.fsp); % call method to generate the PDO.
            end

            %% Project FSP result onto species of interest.
            if obj.useHybrid
                speciesStochastic = setdiff(obj.species,obj.hybridOptions.upstreamODEs);
            else
                speciesStochastic = obj.species;
            end
            Nd = length(speciesStochastic);
            for i=Nd:-1:1
                indsPlots(i) = max(contains(obj.dataSet.linkedSpecies(:,1),speciesStochastic(i)));
            end

            szP = zeros(1,Nd);
            for it = length(obj.tSpan):-1:1
                if ~computeSensitivity||nargout<2
                    szP = max(szP,size(solutions.fsp{it}.p.data));
                else
                    szP = max(szP,size(solutions.sens.data{it}.p.data));
                end
            end

            P = zeros([length(obj.tSpan),szP(indsPlots)]);
            for it = length(obj.tSpan):-1:1
                INDS = setdiff([1:Nd],find(indsPlots));
                if ~computeSensitivity||nargout<2
                    px = solutions.fsp{it}.p;
                else
                    if computeSensitivity&&nargout>=2
                        px = solutions.sens.data{it}.p;
                        Sx = solutions.sens.data{it}.S;
                        parCount = length(Sx);
                        % Add effect of PDO.
                        if ~isempty(obj.pdoOptions.PDO)
                            for iPar = 1:parCount
                                Sx(iPar) = obj.pdoOptions.PDO.computeObservationDistDiff(px, Sx(iPar), iPar);
                            end
                        end
                    end
                end

                % Add effect of PDO.
                if ~isempty(obj.pdoOptions.PDO)
                    try
                        px = obj.pdoOptions.PDO.computeObservationDist(px);
                    catch
                        obj.pdoOptions.PDO = obj.generatePDO(obj.pdoOptions,[],solutions.fsp); % call method to generate the PDO.
                        px = obj.pdoOptions.PDO.computeObservationDist(px);
                    end
                end

                if ~isempty(INDS)
                    d = double(px.sumOver(INDS).data);
                else
                    d = double(px.data);
                end

                P(it,d~=0) = d(d~=0);

                if computeSensitivity&&nargout>=2
                    for iPar = parCount:-1:1
                        if ~isempty(INDS)
                            d = double(Sx(iPar).sumOver(INDS).data);
                        else
                            d = double(Sx(iPar).data);
                        end
                        S{iPar}(it,d~=0) = d(d~=0);
                    end
                end
                %             end
            end

            %% Padd P or Data to match sizes of tensors.
            NP = size(P);
            PD = obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor;
            NDat = size(PD);
            if length(NP)<Nd; NP(end+1:Nd)=1; end
            if max(NDat(2:end)-NP(2:length(NDat)))>0   % Pad if data longer than model
                NP(2:length(NDat)) = max(NP(2:length(NDat)),NDat(2:end));
                tmp = 'P(end';
                for j = 2:length(NDat)
                    tmp = [tmp,',NP(',num2str(j),')'];
                end
                tmp = [tmp,')=0;'];
                eval(tmp);
                if computeSensitivity&&nargout>=2
                    for iPar = 1:parCount
                        tmp2 = strrep(tmp,'P(end',['S{',num2str(iPar),'}(end']);
                        eval(tmp2);
                    end
                end
            end
            if max(NP(2:length(NDat))-NDat(2:end))>0   % Pad if model longer than data
                NDat(2:length(NDat)) = max(NP(2:length(NDat)),NDat(2:end));
                tmp = 'PD(end';
                for j = 2:length(NDat)
                    tmp = [tmp,',NDat(',num2str(j),')'];
                end
                tmp = [tmp,')=0;'];
                eval(tmp);
                %                 if computeSensitivity&&nargout>=2
                %                     for iPar = 1:parCount
                %                         tmp2 = strrep(tmp,'P(end',['S{',num2str(iPar),'}(end']);
                %                         eval(tmp2);
                %                     end
                %                 end
            end
            obj.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor=PD;

            %             if max(NP(2:length(NDat))-NDat(2:end))>0   % truncate if model longer than data
            %                 tmp = 'P = P(:';
            %                 for j = 2:length(NDat)
            %                     tmp = [tmp,',1:',num2str(NDat(j))];
            %                 end
            %                 for j = (length(NDat)+1):4
            %                     tmp = [tmp,',1'];
            %                 end
            %                 tmp = [tmp,');'];
            %                 eval(tmp)
            %                 if computeSensitivity&&nargout>=2
            %                     for iPar = 1:parCount
            %                         tmp2 = strrep(tmp,'P = P',['S{',num2str(iPar),'} = S{',num2str(iPar),'}']);
            %                         eval(tmp2);
            %                     end
            %                 end
            %             end

            P = max(P,1e-10);
            if ~isreal(P)
                P = real(P);
                disp('removed imagionary elements of FSP solution')
            end

            %% Data times for fitting
            if strcmp(obj.fittingOptions.timesToFit,'all')
                times = obj.dataSet.times;
                fitSolutions.ParEstFitTimesList = obj.dataSet.app.ParEstFitTimesList;
                obj.fittingOptions.timesToFit = ones(1,length(obj.dataSet.app.ParEstFitTimesList.Value),'logical');
            else
                times = obj.dataSet.times(obj.fittingOptions.timesToFit);
                fitSolutions.ParEstFitTimesList = obj.dataSet.app.ParEstFitTimesList;
                fitSolutions.ParEstFitTimesList.Value = obj.dataSet.app.ParEstFitTimesList.Value(obj.fittingOptions.timesToFit);
                if ndims(PD)==2
                    PD = PD(find(obj.fittingOptions.timesToFit),:);
                elseif ndims(PD)==3
                    PD = PD(find(obj.fittingOptions.timesToFit),:,:);
                elseif ndims(PD)==4
                    PD = PD(find(obj.fittingOptions.timesToFit),:,:,:);
                end

            end
            %% Compute log likelihood using equal sized P and Data tensors.
            if nargout>=3
                sz = size(P);
                fitSolutions.DataLoadingAndFittingTabOutputs.fitResults.current = zeros([length(times),sz(2:end)]);
                fitSolutions.DataLoadingAndFittingTabOutputs.fitResults.currentData = zeros([length(times),sz(2:end)]);
                fitSolutions.NameTable.Data = [speciesStochastic,speciesStochastic];
                fitSolutions.SpeciesForFitPlot.Value = speciesStochastic(indsPlots);
                fitSolutions.SpeciesForFitPlot.Items = speciesStochastic;
                fitSolutions.DataLoadingAndFittingTabOutputs.dataTensor = PD;
                fitSolutions.FspPrintTimesField.Value = ['[',num2str(obj.tSpan),']'];
                if ~computeSensitivity
                    fitSolutions.FspTabOutputs.solutions = solutions.fsp;
                else
                    fitSolutions.FspTabOutputs.solutions = solutions;
                end
                fitSolutions.FIMTabOutputs.distortionOperator = obj.pdoOptions.PDO;
                fitSolutions.DataLoadingAndFittingTabOutputs.fittingOptions.dataTimes = obj.dataSet.times(obj.fittingOptions.timesToFit);
                fitSolutions.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times = times;
            end

            if computeSensitivity&&nargout>=2
                dlogL_dPar = zeros(parCount,length(times));
            end
            LogLk = zeros(1,length(times));
            KS = zeros(1,length(times));
            numCells = zeros(1,length(times));
            perfectMod = zeros(1,length(times));
            perfectModSmoothed = zeros(1,length(times));
            for i=1:length(times)
                [~,j] = min(abs(obj.tSpan-times(i)));
                if length(times)>1
                    Jind = PD.subs(:,1) == i;
                    SpInds = PD.subs(Jind,:);
                else
                    Jind = ones(size(PD.subs),'logical');
                    SpInds = [ones(length(Jind),1),PD.subs(Jind,:)];
                end
                SpVals = PD.vals(Jind);
                H = sptensor([ones(length(SpVals),1),SpInds(:,2:end)],SpVals,[1,NDat(2:end)]);
                H = double(H);
                Pt = P(j,:,:,:,:,:,:,:,:,:);
                Pt = Pt/max(1,sum(Pt,'all'));
                LogLk(i) = sum(H(:).*log(Pt(:)));
                numCells(i) = sum(H(:));
                if computeSensitivity&&nargout>=2
                    for iPar = parCount:-1:1
                        St = S{iPar}(j,:,:,:,:,:,:);
                        dlogL_dPar(iPar,i) = sum(H(:).*St(:)./Pt(:));
                    end
                end
                if nargout>=3
                    Q = H(:)/sum(H(:));
                    KS(i) = max(abs(cumsum(Q(:))-cumsum(Pt(:))));
                    smQ = smooth(Q);
                    logQ = log(Q); logQ(H==0)=1;
                    logSmQ = log(smQ); logSmQ(H==0)=1;
                    perfectMod(i) = sum(H(:).*logQ);
                    perfectModSmoothed(i) = sum(H(:).*logSmQ);
                    fitSolutions.DataLoadingAndFittingTabOutputs.fitResults.current(i,:,:,:,:,:,:) = Pt;
                    fitSolutions.DataLoadingAndFittingTabOutputs.fitResults.currentData(i,:,:,:,:,:,:) = ...
                        reshape(Q,size(fitSolutions.DataLoadingAndFittingTabOutputs.fitResults.currentData(i,:,:,:)));
                end
            end
            logL = sum(LogLk) + logPrior;
            if imag(logL)~=0
                disp('Imaginary likelihood set to -inf.')
                logL = -inf;
            end
            if nargout>=3
                fitSolutions.DataLoadingAndFittingTabOutputs.V_LogLk = LogLk;
                fitSolutions.DataLoadingAndFittingTabOutputs.numCells = numCells;
                fitSolutions.DataLoadingAndFittingTabOutputs.perfectMod = perfectMod;
                fitSolutions.DataLoadingAndFittingTabOutputs.perfectModSmoothed = perfectModSmoothed;
                fitSolutions.DataLoadingAndFittingTabOutputs.V_KS = KS;
            end
            if computeSensitivity&&nargout>=2
                gradient = sum(dlogL_dPar,2); % need to also add gradient wrt prior!!
            else
                gradient = [];
            end
        end

        function fitErrors = likelihoodSweep(obj,parIndices,scalingRange,makePlot)
            % likelihoodSweep - sweep over range of parameters and return
            % likelihood function values at all parameter combinations.
            arguments
                obj
                parIndices
                scalingRange = linspace(0.5,1.5,15);
                makePlot = false
            end
            obj.fittingOptions.modelVarsToFit = parIndices;  % Choose which parameters to vary.
            pars0 = [obj.parameters{obj.fittingOptions.modelVarsToFit,2}];
            Ngrid=length(scalingRange);
            fitErrors = zeros(Ngrid,Ngrid);
            for i = 1:Ngrid
                for j = 1:Ngrid
                    pars = pars0.*scalingRange([i,j]);
                    fitErrors(i,j) = obj.computeLikelihood(pars);
                end
            end
            if makePlot
                figure
                if length(parIndices)>2
                    disp('plots are only created for first two parameters')
                end
                contourf(scalingRange*pars0(1),scalingRange*pars0(2),fitErrors,30)
                set(gca,'fontsize',15)
                xlabel(obj.parameters{obj.fittingOptions.modelVarsToFit(1)});
                ylabel(obj.parameters{obj.fittingOptions.modelVarsToFit(2)});
                colorbar
                hold on

                [tmp,I] = max(fitErrors);
                [~,J] = max(tmp);
                plot(scalingRange([1,Ngrid])*pars0(1),pars0(2)*[1,1],'k--','linewidth',3)
                plot(pars0(1)*[1,1],scalingRange([1,Ngrid])*pars0(2),'k--','linewidth',3)
                plot(pars0(1)*scalingRange(J),pars0(2)*scalingRange(I(J)),'ro','MarkerSize',20,'MarkerFaceColor','r')
            end
        end

        function [pars,likelihood,otherResults] = maximizeLikelihood(obj,parGuess,fitOptions,fitAlgorithm)
            arguments
                obj
                parGuess = [];
                fitOptions = optimset('Display','iter','MaxIter',10);
                fitAlgorithm = 'fminsearch';
            end

            if isempty(obj.propensitiesGeneral)
                obj = formPropensitiesGeneral(obj);
            end

            if ischar(obj.fittingOptions.modelVarsToFit)&&strcmp(obj.fittingOptions.modelVarsToFit,'all')
                obj.fittingOptions.modelVarsToFit = (1:size(obj.parameters,1));
            end
            if isempty(parGuess)
                parGuess = [obj.parameters{obj.fittingOptions.modelVarsToFit,2}]';
            end

            if strcmp(obj.solutionScheme,'fspSens')   % Set solution scheme to FSP.
                obj.solutionScheme='FSP';
            end

            if strcmp(obj.solutionScheme,'FSP')   % Set solution scheme to FSP.
                [FSPsoln,bounds] = obj.solve;  % Solve the FSP analysis
                obj.fspOptions.bounds = bounds;% Save bound for faster analyses
                objFun = @(x)-obj.computeLikelihood(exp(x),FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.

                if isfield(fitOptions,'suppressFSPExpansion')&&fitOptions.suppressFSPExpansion
                    tmpFSPtol = obj.fspOptions.fspTol;
                    obj.fspOptions.fspTol = inf;
                end
            elseif strcmp(obj.solutionScheme,'ode')  % Set solution scheme to ode.
                objFun = @(x)-obj.computeLikelihoodODE(exp(x));  % We want to MAXIMIZE the likelihood.
            end

            x0 = log(parGuess);

            switch fitAlgorithm
                case 'fminsearch'
                    allFitOptions.suppressFSPExpansion = true;
                    fNames = fieldnames(fitOptions);
                    for i=1:length(fNames)
                        allFitOptions.(fNames{i}) = fitOptions.(fNames{i});
                    end

                    [x0,likelihood,~,otherResults]  = fminsearch(objFun,x0,allFitOptions);

                case 'fminunc'
                    obj.fspOptions.fspTol = inf;
                    objFun = @obj.minusLogL;  % We want to MAXIMIZE the likelihood.
                    x0 = log(parGuess);
                    [x0,likelihood]  = fminunc(objFun,x0,fitOptions,FSPsoln.stateSpace,true);

                case 'particleSwarm'
                    obj.fspOptions.fspTol=inf;
                    rng('shuffle')
                    OBJps = @(x)objFun(x');
                    LB = -5*ones(size(x0'));
                    UB = 5*ones(size(x0'));
                    initSwarm = repmat(x0',fitOptions.SwarmSize-1,1);
                    initSwarm = [x0';initSwarm.*(1+0.1*randn(size(initSwarm)))];
                    fitOptions.InitialSwarmMatrix = initSwarm;
                    [x0,likelihood] = particleswarm(OBJps,length(x0),LB,UB,fitOptions);

                case 'mlSearch'
                    % Not yet working efficiently.
                    allFitOptions.maxIter=1000;
                    allFitOptions.burnIn=30;
                    allFitOptions.updateRate=10;
                    allFitOptions.guessRate=1000;
                    allFitOptions.proposalDistribution=@(x)x+0.01*randn(size(x));
                    allFitOptions.useFIMforSearch = false;
                    allFitOptions.CovFIMscale = 0.6;
                    allFitOptions.suppressFSPExpansion = true;
                    allFitOptions.logForm = true;
                    allFitOptions.plotFunVals = false;
                    allFitOptions.proposalDistributionWide=@(x)x+randn(size(x));

                    fNames = fieldnames(fitOptions);
                    for i=1:length(fNames)
                        allFitOptions.(fNames{i}) = fitOptions.(fNames{i});
                    end

                    if allFitOptions.logForm
                        objFun = @(x)obj.computeLikelihood(exp(x),FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.
                        x0 = log(parGuess);
                    else
                        objFun = @(x)obj.computeLikelihood(x,FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.
                        x0 = (parGuess);
                    end

                    [x0,likelihood]  = mlSearch(objFun,x0,allFitOptions);

                case 'MetropolisHastings'

                    allFitOptions.isPropDistSymmetric=true;
                    allFitOptions.thin=1;
                    allFitOptions.numberOfSamples=1000;
                    allFitOptions.burnIn=100;
                    allFitOptions.progress=true;
                    allFitOptions.proposalDistribution=@(x)x+0.1*randn(size(x));
                    allFitOptions.numChains = 1;
                    allFitOptions.useFIMforMetHast = false;
                    allFitOptions.CovFIMscale = 0.6;
                    allFitOptions.suppressFSPExpansion = true;
                    allFitOptions.logForm = true;

                    j=1;
                    while exist(['TMPmh_',num2str(j),'.mat'],'file')
                        j=j+1;
                    end
                    allFitOptions.saveFile = ['TMPmh_',num2str(j),'.mat'];
                    fNames = fieldnames(fitOptions);
                    for i=1:length(fNames)
                        allFitOptions.(fNames{i}) = fitOptions.(fNames{i});
                    end

                    if allFitOptions.logForm
                        OBJmh = @(x)obj.computeLikelihood(exp(x),FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.
                        x0 = log(parGuess);
                    else
                        OBJmh = @(x)obj.computeLikelihood(x,FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.
                        x0 = (parGuess);
                    end

                    if allFitOptions.useFIMforMetHast
                        TMP = obj;
                        TMP.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
                        [sensSoln] = TMP.solve;  % Solve the sensitivity problem

                        if allFitOptions.logForm
                            fimResults = TMP.computeFIM(sensSoln.sens,'log');
                        else
                            fimResults = TMP.computeFIM(sensSoln.sens);
                        end
                        FIM = TMP.evaluateExperiment(fimResults,TMP.dataSet.nCells);

                        FIMfree = FIM{1}(obj.fittingOptions.modelVarsToFit,obj.fittingOptions.modelVarsToFit);

                        if allFitOptions.logForm&&min(eig(FIMfree))<1
                            disp('Warning -- FIM has one or more small eigenvalues.  Reducing proposal width to 10x in those directions. MH Convergence may be slow.')
                            FIMfree = FIMfree + 1*eye(length(FIMfree));
                        end

                        covFree = FIMfree^-1;
                        covFree = allFitOptions.CovFIMscale*(covFree+covFree')/2;
                        allFitOptions.proposalDistribution=@(x)mvnrnd(x,covFree);
                    end

                    if allFitOptions.suppressFSPExpansion
                        obj.fspOptions.fspTol = inf;
                    end

                    rng('shuffle')
                    if allFitOptions.numChains==1
                        [otherResults.mhSamples,otherResults.mhAcceptance,otherResults.mhValue,x0,likelihood] = ...
                            ssit.parest.metropolisHastingsSample(x0',allFitOptions.numberOfSamples,...
                            'logpdf',OBJmh,'proprnd',allFitOptions.proposalDistribution,...
                            'symmetric',allFitOptions.isPropDistSymmetric,...
                            'thin',allFitOptions.thin,'nchain',1,'burnin',allFitOptions.burnIn,...
                            'progress',allFitOptions.progress,...
                            'saveFileName',allFitOptions.saveFile);
                    else
                        try
                            parpool
                        catch
                        end
                        allFitOptions.progress=0;
                        clear tmpMH*
                        parfor iChain = 1:allFitOptions.numChains
                            [mhSamples, mhAcceptance, mhValue,xbest,fbest] = ...
                                ssit.parest.metropolisHastingsSample(x0',allFitOptions.numberOfSamples,...
                                'logpdf',OBJmh,'proprnd',allFitOptions.proposalDistribution,'symmetric',...
                                allFitOptions.isPropDistSymmetric,...
                                'thin',allFitOptions.thin,'nchain',1,'burnin',allFitOptions.burnIn,...
                                'progress',allFitOptions.progress);
                            tmpMHSamp(iChain) = {mhSamples};
                            tmpMHAcceptance(iChain) = {mhAcceptance};
                            tmpMHValue(iChain) = {mhValue};
                            tmpMHxbest(iChain) = {xbest};
                            tmpMHfbest(iChain) = fbest;
                        end
                        [~,jBest] = max(tmpMHfbest);
                        x0 = tmpMHxbest{jBest}';
                        otherResults.mhSamples = tmpMHSamp;
                        otherResults.mhAcceptance = tmpMHAcceptance;
                        otherResults.mhValue = tmpMHValue;
                        clear tmpMH*
                    end
                    % If fit was in linear space, need to convert to log
                    % space before returning parameters.
                    if ~allFitOptions.logForm
                        pars = log(x0);
                    end

            end

            pars = exp(x0);

            if strcmp(obj.solutionScheme,'FSP')&&isfield(fitOptions,'suppressFSPExpansion')
                obj.fspOptions.fspTol = tmpFSPtol;
            end

        end

        %% Model Reduction Functions
        function [obj,fspSoln] = computeModelReductionTransformMatrices(obj,fspSoln,phi)
            % This function computes linear transformation matrices (PHI
            % and PHIinv) that can be used to switch between the reduced
            % and origional FSP bases.
            arguments
                obj
                fspSoln = []
                phi = []
            end

            numConstraints = length(obj.fspOptions.bounds);

            % Assemble for generator matrix for original FSP problem.
            if ~isfield(fspSoln,'A_total')
                fspSoln.Afsp = ssit.FspMatrix(obj.propensitiesGeneral, [obj.parameters{:,2}]', fspSoln.stateSpace, numConstraints);
                fspSoln.A_total = fspSoln.Afsp.createSingleMatrix(obj.tSpan(1), [obj.parameters{:,2}]');
            end

            % Remove FSP Sinks
            fspSoln.A_total = fspSoln.A_total(1:end-numConstraints,1:end-numConstraints);

            % Call function to compute transformation matrices.
            [obj.modelReductionOptions.phi,...
                obj.modelReductionOptions.phi_inv,...
                obj.modelReductionOptions.phiScale,...
                obj.modelReductionOptions.phiPlot,...
                obj.modelReductionOptions.redOutputs] = ...
                ssit.fsp_model_reduction.getTransformMatrices(...
                obj.modelReductionOptions,...
                fspSoln);

            fspSoln.tOut = obj.tSpan;
            obj.modelReductionOptions.fspSoln=fspSoln;

        end

        function redSolutions = solveReducedFSP(obj)
            arguments
                obj
            end

            numConstraints = length(obj.fspOptions.bounds);
            stateCount = obj.modelReductionOptions.fspSoln.stateSpace.getNumStates();
            % Use Approximate steady state as initial distribution if requested.
            if obj.fspOptions.initApproxSS
                jac = obj.modelReductionOptions.fspSoln.A_total;
                jac = jac(1:end-numConstraints,1:end-numConstraints);
                jac = jac+diag(sum(jac));
                try
                    warning('off')
                    [eigVec,~] = eigs(jac,1,'smallestabs');
                catch
                    try
                        [eigVec,~] = eigs(jac,0);
                    catch
                        try
                            eigVec = null(full(jac));
                        catch
                            disp('Could not find null space. Using uniform.')
                            eigVec = ones(size(jac,1),1);
                        end
                    end
                end
                obj.modelReductionOptions.fspSoln.P0 = [eigVec/sum(eigVec);zeros(numConstraints,1)];
            else % otherwise use user supplied IC.
                obj.modelReductionOptions.fspSoln.P0  = zeros(stateCount + numConstraints, 1);
                obj.modelReductionOptions.fspSoln.P0(1:size(obj.initialCondition,2)) = obj.initialProbs;
            end

            if strcmp(obj.modelReductionOptions.reductionType,'Balanced Model Truncation (HSV)')
                %                 sys = ss(fspSoln.A_total,fspSoln.P0,eye(nStates),[]);
                %                 sysred = balred(sys,n,redOutputs.info);
                %                 A_red = sysred.A;
                %                 q0 = sysred.B;
                %                 OutPutC = sysred.C;
            else
                q0 = obj.modelReductionOptions.phi_inv*obj.modelReductionOptions.fspSoln.P0;
                A_red = obj.modelReductionOptions.phi_inv*...
                    obj.modelReductionOptions.fspSoln.A_total*...
                    obj.modelReductionOptions.phi;
            end

            fspErrorCondition.tInit = obj.modelReductionOptions.fspSoln.tOut(1);
            [~, ~, ~, ~, solutionsNow] = ssit.fsp_ode_solvers.expv_modified(...
                obj.modelReductionOptions.fspSoln.tOut(end), A_red, q0,...
                1e-8, 30,...
                [],...
                obj.modelReductionOptions.fspSoln.tOut,...
                1e-3, [],...
                obj.modelReductionOptions.fspSoln.tOut(1),...
                fspErrorCondition);

            if strcmp(obj.modelReductionOptions.reductionType,'Balanced Model Truncation (HSV)')
                %                 redSolutionsNow = solutionsNow*OutPutC';
            else
                redSolutionsNow = solutionsNow*obj.modelReductionOptions.phi';
                redSolutionsNow = diag(1./sum(redSolutionsNow,2))*redSolutionsNow;
            end

            for j=size(redSolutionsNow,1):-1:1
                redSolutions.fsp{j} = struct(time=obj.modelReductionOptions.fspSoln.tOut(j),...
                    p=ssit.FspVector(obj.modelReductionOptions.fspSoln.stateSpace.states,...
                    redSolutionsNow(j,1:stateCount)),...
                    sinks=[]);
            end
            redSolutions.stateSpace = obj.modelReductionOptions.fspSoln.stateSpace.states;
        end

        %% Plotting/Visualization Functions
        function makePlot(obj,solution,plotType,indTimes,includePDO,figureNums,lineProps)
            % SSIT.makePlot -- tool to make plot of the FSP or SSA results.
            % arguments:
            %   solution -- solution structure from SSIT.
            %   plotType - chosen type of plot:
            %       FSP options: 'means' -- mean versus time
            %                    'meansAndDevs' -- means +/- STD vs time
            %                    'marginals' -- marginal distributions over
            %                           time
            %                    'joints' -- joint distributions vs time.
            %       sensFSP options:
            %                   'marginals' -- sensitivity of marginal distributions
            %                           for each parameter and time point.
            %       SSA options: 'means' -- mean versus time
            %                    'meansAndDevs' -- means +/- STD vs time
            %                    'trajectories' -- set of individual trajectories vs time.
            %
            % examples:
            %
            %   F = SSIT('ToggleSwitch')
            %   F.solutionScheme = 'FSP'
            %   [FSPsoln,bounds] = F.solve;  % Returns the solution and the
            %                             % bounds for the FSP projection
            %   F.makePlot(FSPsoln,'marginals')  % Make plot of FSP
            %                                    % marginal distributions.

            %   F.solutionScheme = 'fspSens'
            %   [sensSoln,bounds] = F.solve;  % Returns the sensitivity and the
            %                                 bounds for the FSP projection
            %   F.makePlot(sensSoln,'marginals')% Make plot of
            %                                   %sensitivities of marginal
            %                                   distributions at final
            %                                   time.
            arguments
                obj
                solution
                plotType = 'means';
                indTimes = [];
                includePDO = false;
                figureNums = [];
                lineProps = {'linewidth',2};
            end
            if isempty(figureNums)
                h =  findobj('type','figure');
                if isfield(h,'Number')
                    figureNums = max([h.Number])+(1:10);
                else
                    figureNums = (1:10);
                end      
            end
            kfig = 1;
            switch obj.solutionScheme
                case 'FSP'
                    app.FspTabOutputs.solutions = solution.fsp;
                    if includePDO
                        if ~isempty(obj.pdoOptions.PDO)
                            for i=1:length(app.FspTabOutputs.solutions)
                                app.FspTabOutputs.solutions{i}.p = obj.pdoOptions.PDO.computeObservationDist(app.FspTabOutputs.solutions{i}.p);
                            end
                        else
                            warning('obj.pdoOptions.PDO has not been set')
                        end
                    end
                    app.FspPrintTimesField.Value = mat2str(obj.tSpan);
                    solution = exportFSPResults(app);
                    Nd = length(solution.Marginals{end});
                    if isempty(indTimes)
                        indTimes = 1:length(solution.T_array);
                    end
                    Nt = length(indTimes);
                    switch plotType
                        case 'means'
                            plot(solution.T_array(indTimes),solution.Means(indTimes,:),lineProps{:});
                        case 'meansAndDevs'
                            figure(figureNums(kfig)); kfig=kfig+1;
                            for i = 1:Nd
                                subplot(Nd,1,i); hold on
                                errorbar(solution.T_array(indTimes),solution.Means(indTimes,i),sqrt(solution.Var(indTimes,i)),lineProps{:});
                                ylabel(obj.species{i})
                            end
                            xlabel('time')
                        case 'marginals'
                            for j = 1:Nd
                                f = figure(figureNums(kfig)); kfig=kfig+1;
                                f.Name = ['Marginal Distributions of ',obj.species{j}];
                                Nr = ceil(sqrt(Nt));
                                Nc = ceil(Nt/Nr);
                                for i = 1:Nt
                                    i2 = indTimes(i);
                                    subplot(Nr,Nc,i); hold on
                                    stairs(solution.Marginals{i2}{j},lineProps{:});
                                    set(gca,'fontsize',15)
                                    title(['t = ',num2str(solution.T_array(i2),2)])
                                end
                            end
                        case 'joints'
                            if Nd<2
                                error('Joint distributions only avaialble for 2 or more species.')
                            else
                                for j1 = 1:Nd
                                    for j2 = j1+1:Nd
                                        h = figure(figureNums(kfig)); kfig=kfig+1;
                                        h.Name = ['Joint Distribution of ',obj.species{j1},' and ',obj.species{j2}];
                                        Nr = ceil(sqrt(Nt));
                                        Nc = ceil(Nt/Nr);
                                        for i = 1:Nt
                                            i2 = indTimes(i);
                                            subplot(Nr,Nc,i);
                                            contourf(log10(max(1e-16,solution.Joints{i2}{j1,j2})));
                                            colorbar
                                            title(['t = ',num2str(solution.T_array(i2),2)])
                                            %                                             if mod(i-1,Nc)==0;
                                            ylabel(['x',num2str(j1)]);
                                            %                                             end
                                            %                                             if (i+Nc)>Nt;
                                            xlabel(['x',num2str(j2)]);
                                            %                                             end
                                            set(gca,'FontSize',15)
                                        end
                                    end
                                end
                            end
                        case 'escapeTimes'
                            f = figure(figureNums(kfig)); kfig=kfig+1;
                            subplot(2,1,1)
                            z = solution.EscapeCDF(indTimes,:);
                            t = solution.T_array(indTimes);
                            plot(t,z,'linewidth',3); hold on
                            set(gca,'fontsize',16)
                            ylabel('Escape CDF')

                            subplot(2,1,2)
                            zp = (z(2:end,:)-z(1:end-1,:))./(t(2:end)-t(1:end-1))';
                            tp = (t(2:end)+t(1:end-1))/2;
                            plot(tp',zp,'linewidth',3); hold on
                            set(gca,'fontsize',16)
                            ylabel('Escape PDF')
                            xlabel('time')
                    end
                case 'SSA'
                    Nd = size(solution.trajs,1);
                    if isempty(indTimes)
                        indTimes = 1:length(solution.T_array);
                    end
                    switch plotType
                        case 'trajectories'
                            figure(figureNums(kfig)); kfig=kfig+1;
                            for i=1:Nd
                                subplot(Nd,1,i)
                                plot(solution.T_array(indTimes),squeeze(solution.trajs(i,indTimes,:)));
                            end
                        case 'means'
                            figure(figureNums(kfig)); kfig=kfig+1;
                            plot(solution.T_array(indTimes),squeeze(mean(solution.trajs(:,indTimes,:),3)));
                        case 'meansAndDevs'
                            figure(figureNums(kfig)); kfig=kfig+1;
                            vars = var(solution.trajs(:,indTimes,:),[],3);
                            errorbar(solution.T_array(indTimes),squeeze(mean(solution.trajs(:,indTimes,:),3)),sqrt(vars));
                    end
                case 'fspSens'

                    if includePDO
                        if ~isempty(obj.pdoOptions.PDO)
                            for i=1:length(solution.sens.data)
                                for j=1:length(solution.sens.data{i}.S)
                                    solution.sens.data{i}.S(j) = obj.pdoOptions.PDO.computeObservationDist(solution.sens.data{i}.S(j));
                                end
                            end
                        else
                            warning('obj.pdoOptions.PDO has not been set')
                        end
                    end

                    app.SensFspTabOutputs.solutions = solution.sens;
                    app.SensPrintTimesEditField.Value = mat2str(obj.tSpan);
                    if ~isempty(obj.parameters)
                        app.ReactionsTabOutputs.parameters = obj.parameters(:,1);
                    else
                        app.ReactionsTabOutputs.parameters = [];
                    end
                    app.ReactionsTabOutputs.varNames = obj.species;
                    solution.plotable = exportSensResults(app);

                    Np = size(solution.plotable.sensmdist,1);
                    Nd = size(solution.plotable.sensmdist,2);
                    if isempty(indTimes)
                        indTimes = length(solution.plotable.T_array);
                    end
                    Nt = length(indTimes);
                    Nr = ceil(sqrt(Np));
                    Nc = ceil(Np/Nr);
                    switch plotType
                        case 'marginals'
                            for it = 1:Nt
                                it2 = indTimes(it);
                                for id = 1:Nd
                                    f = figure(figureNums(kfig)); kfig=kfig+1;
                                    f.Name = ['Marg. Dist. Sensitivities of x',num2str(id),' at t=',num2str(solution.plotable.T_array(it2))];
                                    for j = 1:Np
                                        subplot(Nr,Nc,j); hold on;
                                        stairs(solution.plotable.sensmdist{j,id,it2},lineProps{:});
                                        set(gca,'fontsize',15)
                                        title(obj.parameters{j,1})
                                        %                                         if mod(j-1,Nc)==0;
                                        ylabel(['sensitivity']);
                                        %                                         end
                                        %                                         if (j+Nc)>Np;
                                        xlabel(['x',num2str(id)]);
                                        %                                         end
                                    end
                                end
                            end
                    end
            end
        end

        function makeFitPlot(obj,fitSolution,smoothWindow,fignums)
            % Produces plots to compare model to experimental data.
            arguments
                obj
                fitSolution =[];
                smoothWindow = 5;
                fignums = [];
            end
            if isempty(fitSolution)
                [~,~,fitSolution] = obj.computeLikelihood;
            end
            makeSeparatePlotOfData(fitSolution,smoothWindow,fignums)
        end

        function plotMHResults(obj,mhResults,FIM)
            arguments
                obj
                mhResults = [];
                FIM =[];
            end

            parNames = obj.parameters(obj.fittingOptions.modelVarsToFit,1);
            if ~isempty(FIM)
                pars = [obj.parameters{obj.fittingOptions.modelVarsToFit,2}];
               
                parsLog = log10(pars);

                if ~iscell(FIM)
                    FIM = diag(pars)*...
                        FIM(obj.fittingOptions.modelVarsToFit,obj.fittingOptions.modelVarsToFit)*...
                        diag(pars);
                    covFIM{1} = FIM^(-1)/log(10)^2;
                else
                    for i=1:length(FIM)
                        FIMi = diag(pars)*...
                            FIM{i}(obj.fittingOptions.modelVarsToFit,obj.fittingOptions.modelVarsToFit)*...
                            diag(pars);
                        covFIM{i} = FIMi^(-1)/log(10)^2;
                    end
                end
            end


            if ~isempty(mhResults)
                % Make figures for MH convergence
                figure;
                plot(mhResults.mhValue);
                xlabel('Iteration number');
                ylabel('log-likelihood')
                title('MH Convergence')

                figure
                ac = xcorr(mhResults.mhValue-mean(mhResults.mhValue),'normalized');
                ac = ac(size(mhResults.mhValue,1):end);
                plot(ac,'LineWidth',3); hold on
                N = size(mhResults.mhValue,1);
                tau = 1+2*sum((ac(2:N/100)));
                Neff = N/tau;
                xlabel('Lag');
                ylabel('Auto-correlation')
                title('MH Convergence')

                figure
                [valDoneSorted,J] = sort(mhResults.mhValue);
                smplDone = mhResults.mhSamples(J,:);
                Np = size(mhResults.mhSamples,2);
            end
            
            fimCols = {'k','c','b','g','r'};

            for i=1:Np-1
                for j = i+1:Np
                    subplot(Np-1,Np-1,(i-1)*(Np-1)+j-1);

                    if ~isempty(mhResults)
                        scatter(smplDone(:,j)/log(10),smplDone(:,i)/log(10),20,valDoneSorted,'filled'); hold on;
                        par0 = mean(smplDone(:,[j,i])/log(10));
                        cov12 = cov(smplDone(:,j)/log(10),smplDone(:,i)/log(10));
                    end
                    if ~isempty(FIM)
                        for iFIM = 1:length(covFIM)
                            ssit.parest.ellipse(parsLog([j,i]),icdf('chi2',0.9,2)*covFIM{iFIM}([j,i],[j,i]),fimCols{mod(iFIM,5)+1},'linewidth',2)
                        end
                    end
                    if ~isempty(mhResults)
                        ssit.parest.ellipse(par0,icdf('chi2',0.9,2)*cov12,'m--','linewidth',2)
                    end
                    xlabel(['log_{10}(',parNames{j},')']);
                    ylabel(parNames{i});
                end
            end
        end

        function makeMleFimPlot(obj,MLE,FIM,indPars,CI,figNum,par0)
            arguments
                obj
                MLE = []
                FIM = []
                indPars = [1,2];
                CI = 0.95
                figNum=[]
                par0 = []
            end
            if isempty(figNum)
                gcf;
            end

            CIp = round(CI*100);

            legs = {};

            if ~isempty(MLE)
                scatter(MLE(indPars(1),:),MLE(indPars(2),:),100*ones(size(MLE(indPars(1),:))),'filled');
                covMLE = cov(MLE');
                muMLE = mean(MLE,2);
                hold on
                ssit.parest.ellipse(muMLE(indPars),icdf('chi2',CI,2)*covMLE(indPars,indPars),'linewidth',3)
                legs(end+1:end+2) = {['MLE, N=',num2str(length(MLE))],[num2str(CIp),'% CI (MLE)']};
                if isempty(par0)
                    par0 = muMLE;
                end
            end

            if ~isempty(FIM)
                covFIM = FIM^(-1);
                ssit.parest.ellipse(par0(indPars),icdf('chi2',CI,2)*covFIM(indPars,indPars),'--','linewidth',3)
                legs(end+1) = {[num2str(CIp),'% CI (FIM)']};
            end
            set(gca,'fontsize',15)
            legend(legs)

        end

    end
    methods (Static)
        function FIM = totalFim(fims,Nc)
            Nt = size(fims,1);
            Ns = size(fims,2);
            FIM = cell(1,Ns);
            for is = 1:Ns
                FIM{is} = 0*fims{1};
                for it = 1:Nt
                    FIM{is} = FIM{is}+Nc(it)*fims{it,is};
                end
            end
        end
        function k = findBestMove(fims,Ncp,met,NcMax)
            arguments
                fims
                Ncp
                met
                NcMax = [];
            end
            Nt = size(fims,1);
            if isempty(NcMax)
                NcMax = inf*ones(1,Nt);
            end
            Ns = size(fims,2);
            obj = zeros(Nt,Ns);
            FIM0 = SSIT.totalFim(fims,Ncp);
            for is = 1:Ns
                for it = 1:Nt
                    if Ncp(it)<NcMax(it)
                        % If one can do that experiment.
                        FIM = FIM0{is}+fims{it,is};
                    else
                        % If there are no more cells avalable for that time
                        % point.
                        FIM = FIM0{is};
                    end
                    obj(it,is) = met(FIM);
                end
            end
            [~,k] = min(mean(obj,2));
        end
        function propensities = parameterizePropensities(GenProps,Parset)
            propensities = GenProps;
            for i=1:length(GenProps)
                if ~isempty(propensities{i}.stateDependentFactor)
                    propensities{i}.stateDependentFactor = @(x)GenProps{i}.stateDependentFactor(x,Parset);
                end
                if ~isempty(propensities{i}.hybridFactor)
                    propensities{i}.hybridFactor = @(t,v)GenProps{i}.hybridFactor(t,Parset,v');
                end
                if ~isempty(propensities{i}.hybridFactorVector)
                    propensities{i}.hybridFactorVector = @(t,v)GenProps{i}.hybridFactorVector(t,Parset,v');
                end
                if ~isempty(propensities{i}.timeDependentFactor)
                    propensities{i}.timeDependentFactor = @(t)GenProps{i}.timeDependentFactor(t,Parset);
                end
                if ~isempty(propensities{i}.hybridJointFactor)
                    propensities{i}.hybridJointFactor = @(t,x,v)GenProps{i}.hybridJointFactor(t,x,Parset,v');
                end
            end
        end
    end
end