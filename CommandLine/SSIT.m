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
            'initApproxSS',false); % Options for FSP solver.
        sensOptions = struct('solutionMethod','forward'); 
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
        dataSet = [];
    end

    properties (Dependent)
        fspConstraints % FSP Constraint Functions
        pars_container % Container for parameters
        propensities % Processed propensity functions for use in solvers
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
            %SSIT Construct an instance of the SSIT class
            addpath(genpath('../src'));
            if isempty(modelFile)
                return
            else
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

        function Propensities = get.propensities(obj)
            switch obj.solutionScheme
                case 'FSP'
                    propenType = "str";
                case 'SSA'
                    propenType = "fun";
            end
            propstrings = ssit.SrnModel.processPropensityStrings(obj.propensityFunctions,...
                obj.inputExpressions,...
                obj.pars_container,...
                propenType,...
                obj.species);

            if strcmp(obj.solutionScheme,'FSP')
                n_reactions = length(obj.propensityFunctions);
                Propensities = cell(n_reactions, 1);
                for i = 1:n_reactions
                    Propensities{i} = ssit.Propensity.createFromString(propstrings{i}, obj.stoichiometry(:,i), i);
                end
            elseif strcmp(obj.solutionScheme,'SSA')
                Propensities = propstrings;
            end
        end

        function constraints = get.fspConstraints(obj)
            % Makes a list of FSP constraints that can be used by the FSP
            % solver.
            nSpecies = length(obj.species);
            Data = cell(nSpecies*2,3);
            for i = 1:nSpecies
%                 Data(i,:) = {['-',obj.species{i}],'<',0};
%                 Data(nSpecies+i,:) = {obj.species{i},'<',1};
                Data(i,:) = {['-x',num2str(i)],'<',0};
                Data(nSpecies+i,:) = {['x',num2str(i)],'<',1};
            end
            for i = 1:length(obj.customConstraintFuns)
                Data(2*nSpecies+i,:) = {obj.customConstraintFuns{i},'<',1};
            end
            constraints.f = readConstraintsForAdaptiveFsp([], obj.species, Data);
            if isempty(obj.fspOptions.bounds)
                constraints.b = [Data{:,3}]';
            else
                constraints.b = obj.fspOptions.bounds;
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
            options = optimset('display','iter');
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
        
        %% Model Analysis Functions
        function [Solution, bConstraints] = solve(obj,stateSpace,saveFile)
            arguments
                obj
                stateSpace = [];
                saveFile=[];
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

            switch obj.solutionScheme
                case 'FSP'
                    if ~isempty(stateSpace)&&size(stateSpace.states,2)~=stateSpace.state2indMap.Count
                        error('HERE')
                    end

                    [Solution.fsp, bConstraints,Solution.stateSpace] = ssit.fsp.adaptiveFspSolve(obj.tSpan,...
                        obj.initialCondition,...
                        obj.initialProbs,...
                        obj.stoichiometry, ...
                        obj.propensities, ...
                        obj.fspOptions.fspTol, ...
                        obj.fspConstraints.f, ...
                        obj.fspConstraints.b,...
                        obj.fspOptions.verbose, ...
                        obj.fspOptions.fspIntegratorRelTol, ...
                        obj.fspOptions.fspIntegratorAbsTol, ...
                        obj.fspOptions.odeSolver,stateSpace,...
                        obj.fspOptions.usePiecewiseFSP,...
                        obj.fspOptions.initApproxSS,...
                        obj.species);
                case 'SSA'
                    Solution.T_array = obj.tSpan;
                    Nt = length(Solution.T_array);
                    nSims = obj.ssaOptions.Nexp*obj.ssaOptions.nSimsPerExpt*Nt;
                    W = obj.propensities;
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
                                obj.ssaOptions.signalUpdateRate)
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
                                obj.ssaOptions.signalUpdateRate);
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
                        obj.species);
                    %                     app.SensFspTabOutputs.solutions = Solution.sens;
                    %                     app.SensPrintTimesEditField.Value = mat2str(obj.tSpan);
                    %                     Solution.plotable = exportSensResults(app);
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
                                     A.(['exp',num2str(ie),'_s',num2str(s),'_Distorted'])((it-1)*obj.ssaOptions.nSimsPerExpt+1:obj.ssaOptions.nSimsPerExpt) = ...
                                         Solution.trajsDistorted(s,it,(ie-1)*obj.ssaOptions.nSimsPerExpt+1:obj.ssaOptions.nSimsPerExpt);
                                 end
                             end
%                          end
                     end
                 end
                 writetable(A,saveFile)
                 disp(['FSP Samples saved to ',saveFile])
             end
        end

        function [fimResults,sensSoln] = computeFIM(obj,sensSoln)
            % computeFIM - computes FIM at all time points.
            % Arguments:
            %   sensSoln - (optional) previously compute FSP Sensitivity.
            %              Automatically computed if not provided.
            % Outputs:
            %   fimResults - FIM at each time point in obj.tSpan
            %   sensSoln - FSP Sensitivity.
            arguments
                obj
                sensSoln = [];
            end
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
        end

        function [fimTotal,mleCovEstimate,fimMetrics] = evaluateExperiment(obj,...
                fimResults,cellCounts)
            fimTotal = 0*fimResults{1};
            Np = size(fimTotal,1);
            for i=1:length(cellCounts)
                fimTotal = fimTotal + cellCounts(i)*fimResults{i};
            end

            if nargout>=2
                % Estimate MLE covariance
                if rank(fimTotal)<Np
                    disp(['FIM has rank ',num2str(rank(fimTotal)),' and is not invertable for this experiment design'])
                    mleCovEstimate = NaN;
                else
                    mleCovEstimate = fimTotal^-1;
                end
            end

            if nargout>=3
                % Compute FIM metrics.
                fimMetrics.det = det(fimTotal);
                fimMetrics.trace = trace(fimTotal);
                fimMetrics.minEigVal = min(eig(fimTotal));
            end
        end

        function [Nc] = optimizeCellCounts(obj,fims,nCellsTotal,FIMMetric,Nc)
           % This function optimizes the number of cells per time point
           % according to the user-provide metric.  
           % Allowable metric are:
           %   'Determinant' - maximize the determinant of the FIM
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
            arguments
                obj
                fims
                nCellsTotal
                FIMMetric = 'Smallest Eigenvalue';
                Nc = [];
            end
            switch FIMMetric
                case 'Determinant'
                    met = @(A)-det(A);
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
            NT = size(fims(:,1),1);
            
            if isempty(Nc)
                Nc = zeros(1,NT);
                Nc(1)=nCellsTotal;
            end

            Converged = 0;
            while Converged==0
                Converged = 1;
                for i = 1:NT
                    while Nc(i)>0
                        Ncp = Nc;
                        Ncp(i) = Ncp(i)-1;
                        k = SSIT.findBestMove(fims(:,1),Ncp,met);
                        if k==i
                            break
                        end
                        Nc = Ncp;
                        Nc(k)=Nc(k)+1;
                        Converged = 0;
                    end
                end
            end
        end

        %% Data Loading and Fitting
        function [obj] = loadData(obj,dataFileName,linkedSpecies,conditions)
            arguments
                obj
                dataFileName
                linkedSpecies
                conditions = {};
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
            obj.dataSet.times = unique([obj.dataSet.DATA{:,Jtime}]);

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
            obj.dataSet.app = filterAndMarginalize([],[],obj.dataSet.app);

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

        end

        function [logL,gradient] = minusLogL(obj,pars,stateSpace,computeSensitivity)
            [logL,gradient] = computeLikelihood(obj,exp(pars),stateSpace,computeSensitivity);
            logL = -logL;
            gradient = -gradient.*exp(pars);
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
%                 try
                    [solutions] = obj.solve(stateSpace);  % Solve the FSP analysis
%                 catch
                    % sometimes the FSP analysis crashes if the provided statespace
                    % is incorrect.  Not sure why, but regenerating it will fix
                    % the problem...
                    %                 warning('Regenerate stateset.')
%                     [solutions,bounds] = obj.solve;  % Solve the FSP analysis
%                 end
            end
            obj.parameters =  originalPars;

            if ~isempty(pars)&&nPdoPars>0
                obj.pdoOptions.props.ParameterGuess(indsPdoParsToFit) = pars(nModelPars+1:end);
                obj.pdoOptions.PDO = obj.generatePDO(obj.pdoOptions,[],solutions.fsp); % call method to generate the PDO.
            end

            %% Project FSP result onto species of interest.
            Nd = length(obj.species);
            for i=Nd:-1:1
                indsPlots(i) = max(contains(obj.dataSet.linkedSpecies(:,1),obj.species(i)));
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
                %                 if ~isempty(solutions.fsp{it})
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
                fitSolutions.NameTable.Data = [obj.species,obj.species];
                fitSolutions.SpeciesForFitPlot.Value = obj.species(indsPlots);
                fitSolutions.SpeciesForFitPlot.Items = obj.species;
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
            %             obj.parameters = originalPars;
            %             fit_error = app.DataLoadingAndFittingTabOutputs.J_LogLk;
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
            if isempty(parGuess)
                parGuess = [obj.parameters{obj.fittingOptions.modelVarsToFit,2}]';
            end
            obj.solutionScheme = 'FSP';   % Set solution scheme to FSP.
            [FSPsoln,bounds] = obj.solve;  % Solve the FSP analysis
            obj.fspOptions.bounds = bounds;% Save bound for faster analyses
            objFun = @(x)-obj.computeLikelihood(exp(x),FSPsoln.stateSpace);  % We want to MAXIMIZE the likelihood.
            x0 = log(parGuess);

            switch fitAlgorithm
                case 'fminsearch'
                    [x0,likelihood,~,otherResults]  = fminsearch(objFun,x0,fitOptions);

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
                    allFitOptions.useFIMforMH = false;
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

                        fimResults = TMP.computeFIM(sensSoln.sens);
                        [FIM] = TMP.evaluateExperiment(fimResults,TMP.dataSet.nCells);

                        if allFitOptions.logForm
                            FIMlog = diag([TMP.parameters{obj.fittingOptions.modelVarsToFit,2}])*...
                                FIM(obj.fittingOptions.modelVarsToFit,obj.fittingOptions.modelVarsToFit)*...
                                diag([TMP.parameters{obj.fittingOptions.modelVarsToFit,2}]);
                        else
                            FIMlog = FIM(obj.fittingOptions.modelVarsToFit,obj.fittingOptions.modelVarsToFit);
                        end
                        
                        covLog = FIMlog^-1;
                        covLog = allFitOptions.CovFIMscale*(covLog+covLog')/2;
                        allFitOptions.proposalDistribution=@(x)mvnrnd(x,covLog);                       
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
        end
        %% Plotting/Visualization Functions
        function makePlot(obj,solution,plotType,indTimes,includePDO,figureNums)
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
                figureNums = [1:100];
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
                            plot(solution.T_array(indTimes),solution.Means(indTimes,:));
                        case 'meansAndDevs'
                            figure(figureNums(kfig)); kfig=kfig+1;
                            for i = 1:Nd
                                subplot(Nd,1,i); hold on
                                errorbar(solution.T_array(indTimes),solution.Means(indTimes,i),sqrt(solution.Var(indTimes,i)),'linewidth',2);
                            end
                        case 'marginals'
                            for j = 1:Nd
                                f = figure(figureNums(kfig)); kfig=kfig+1;
                                f.Name = ['Marginal Distributions of x',num2str(j)];
                                Nr = ceil(sqrt(Nt));
                                Nc = ceil(Nt/Nr);
                                for i = 1:Nt
                                    i2 = indTimes(i);
                                    subplot(Nr,Nc,i); hold on
                                    stairs(solution.Marginals{i2}{j},'linewidth',2);
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
                                        h.Name = ['Joint Distribution of x',num2str(j1),' and x',num2str(j2)];
                                        Nr = ceil(sqrt(Nt));
                                        Nc = ceil(Nt/Nr);
                                        for i = 1:Nt
                                            i2 = indTimes(i);
                                            subplot(Nr,Nc,i);
                                            contourf(log10(solution.Joints{i2}{j1,j2}));
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
                                        stairs(solution.plotable.sensmdist{j,id,it2},'linewidth',2);
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
                mhResults 
                FIM =[];
            end

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


            figure
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

            fimCols = {'k','c','b','g','r'};
            
            for i=1:Np-1
                for j = i+1:Np
                    subplot(Np-1,Np-1,(i-1)*(Np-1)+j-1);
                    scatter(smplDone(:,j)/log(10),smplDone(:,i)/log(10),20,valDoneSorted,'filled'); hold on;
                    if ~isempty(FIM)
                        for iFIM = 1:length(covFIM)
                            ssit.parest.ellipse(parsLog([j,i]),icdf('chi2',0.9,2)*covFIM{iFIM}([j,i],[j,i]),fimCols{iFIM},'linewidth',2)
                        end
                    end
                    par0 = mean(smplDone(:,[j,i])/log(10));
                    cov12 = cov(smplDone(:,j)/log(10),smplDone(:,i)/log(10));
                    ssit.parest.ellipse(par0,icdf('chi2',0.9,2)*cov12,'m--','linewidth',2)
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
            FIM = 0*fims{1};
            for i = 1:length(fims)
                FIM = FIM+Nc(i)*fims{i};
            end
        end
        function k = findBestMove(fims,Ncp,met)
            obj = zeros(1,length(Ncp));
            FIM0 = SSIT.totalFim(fims,Ncp);
            for i = 1:length(Ncp)
                FIM = FIM0+fims{i};
                obj(i) = met(FIM);
            end
            [~,k] = min(obj);
        end
    end
end