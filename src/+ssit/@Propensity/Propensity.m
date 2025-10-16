classdef Propensity
    % Object for storing the propensity function of a reaction
    % in stochastic chemical reaction network model input from SSIT.
    %
    % Parameters
    % ----------
    %
    %   reactionIndex:  int
    %       index of the reaction associated with this propensity.
    %
    %   isTimeDependent: logical
    %       true iff this propensity is time-dependent.
    %
    %   isFactorizable: logical
    %       true if this propensity is factorizable into a time-dependent
    %       factor and a state-dependent factor.
    %
    %   stoichVector: column vector
    %       the stoichiometry vector associated with the reaction.
    %
    %   timeDependentFactor: function handle
    %       if isFactorizable == true, the time-dependent factor of the
    %       propensity. This function object must be callable with the
    %       syntax
    %           ``y = f(t)``
    %
    %   stateDependentFactor: function handle
    %       if isFactorizable == true, the state-dependent factor of the
    %       propensity. This function object must be able to take as input
    %       an array of states. In particular, if the input is an array X,
    %       the output must be a row vector Y with length(Y) == size(X, 2).
    %
    %   jointDependentFactor: function handle
    %       if isFactorizable == false, this property must be set. It must
    %       be callable with the syntax
    %           ``Y = f(t, X)``
    %       where t is a scalar representing the time variable, and X a 2-D
    %       array representing states as its columns. The output Y must be
    %       a row vector with length(Y) == size(X, 2).
    %
    %   hybridFactor: function handle
    %       if isFactorizable == true, the hybrid factor of the
    %       propensity. This function object must be callable with the
    %       syntax
    %           ``y = f(t,v)`` where 'v' is the vector of the ODE species.
    %
    %   jointHybridFactor: function handle
    %       if isFactorizable == false and useHybrid == true, this property
    %       must be set. It must be callable with the syntax
    %           ``Y = f(t, X, v)``
    %       where t is a scalar representing the time variable, and X a 2-D
    %       array representing states as its columns, and v is the vector
    %       of the upstream ODE species. The output Y must be
    %       a row vector with length(Y) == size(X, 2).
    properties
        reactionIndex;
        isTimeDependent; % (true/false) telling us if this propensity is time-dep
        isFactorizable; % (true/false) if time-dep, is it isFactorizable?
        stoichVector; % (column vec) the associated stoichiometry vector
        
        timeDependentFactor; % function handle to compute the time-dependent factor,
                % left empty if it is either independent of time or not isFactorizable
        stateDependentFactor; % function handle to compute the state-dependent factor.
        jointDependentFactor; % function handle to compute the time-and-state-dependent factor, this
                % is only used if the propensity is not isFactorizable.
        
        hybridFactor % function handle for factorized hybrid propensity functions.
        hybridFactorVector  % function handle for all hybrid or t-dependent reactions
        xFactorVector % function handle for all state-dependent reactions
        hybridJointFactor % function handle for joint hybrid propensity functions.

        sensTimeFactorVec
        sensStateFactor
        sensStateFactorVec
        
        originalString
       
        parameters % list of required model parameters
        ODEstoichVector=[]; % (column vec) the associated stoichiometry vector (for upstream ODEs)
       
        anonymousT = false; % are time dependent functions anonymous?
        anonymousX = true; % are state dependent functions anonymous?
        DstateFactorDstate % function handle for jacobian with respect to the states.
        DhybridFactorDodesVec % function handle for jacobian with respect to the upstream odes.

    end

    methods (Access = public)
        function obj = Propensity(stoichVector, reactionIndex)
            obj.reactionIndex = reactionIndex;
            obj.stoichVector = stoichVector;
            obj.isTimeDependent = false;
            obj.isFactorizable = true;
        end

        function y = evaluate(obj, t, x, parameters, upstreamSpecies)
            arguments
                obj
                t
                x
                parameters
                upstreamSpecies = [];
            end
            if (size(x, 1) ~= size(obj.stoichVector, 1))
                error('Propensity.eval : input state dimension is not consistent with the stoichiometry vector in the object''s record');
            end
            if (~obj.isTimeDependent)
                y = obj.stateDependentFactor(x);
            else
                if (obj.isFactorizable)
                    y = obj.timeDependentFactor(t)*obj.stateDependentFactor(x,parameters);
                else
                    if ~isempty(obj.jointDependentFactor)
                        y = obj.jointDependentFactor(t, x, parameters);
                    elseif ~isempty(obj.hybridJointFactor)
                        y = obj.hybridJointFactor(t, x, parameters, upstreamSpecies);
                    end
                end
            end
        end

        function y = evaluateStateFactor(obj, x, parameters, varNames, computeSens, ipar)
            arguments
                obj
                x
                parameters
                varNames=[];
                computeSens = false
                ipar  = 1;
            end
            if isempty(varNames)
                varNames={'x1','x2','x3','x4'};
            end
            if (size(x, 1) ~= size(obj.stoichVector, 1))
                error('Propensity.evalAt : input state dimension is not consistent with the stoichiometry vector in the object''s record');
            end
            if computeSens
                if ~isempty(parameters)
                    y = obj.sensStateFactor{ipar}(x,parameters)+zeros(length(parameters),size(x,2));
                    y = y(ipar,:);
                else
                    y = zeros(length(parameters),size(x,2));
                end
            else
                if ~isempty(parameters)
                    y = obj.stateDependentFactor(x,parameters);
                else
                    y = obj.stateDependentFactor(x);
                end
            end
        end
    end
    methods (Access = public, Static)

        function obj = createAsHybridVec(symbolicExpression, stoichMatrix, nonXTpars, species, upstreamODEs, logicTerms, prefixName, computeSens)
            arguments
                symbolicExpression
                stoichMatrix
                nonXTpars
                species
                upstreamODEs = {};
                logicTerms = {};
                prefixName = 'default';
                computeSens = false;
            end

            % Construct an vector of Hybrid Propensity from a symbolic expression.
            %
            % Parameters
            % ----------
            %
            %   symbolicExpression: sym
            %       is a symbolic expression of the propensity.
            %
            %   stoichVector: column vector
            %       is the column stoichiometry vector associated
            %       with this propensity.
            %
            %   reactionIndex: integer
            %       index of the reaction associated with this propensity.
            %
            %   nonXTpars: cell of strings
            %       names of the model parameters (excluding the species
            %       names)
            %
            %   species: cell of strings
            %       names of the species.
            %
            %   upstreamODEs: cell of strings
            %       names of the upstream species that are to be treated as
            %       ODEs.
            %
            % Returns
            % -------
            %
            %   obj: Propensity
            %       an instance of the Propensity class.
            %

            if ~isempty(upstreamODEs)
                jStochastic = find(~contains(species,upstreamODEs));
                jODE = find(contains(species,upstreamODEs));
                timeFunName = 'hybridFactor';
                jntFactorName = 'hybridJointFactor';
            else
                jStochastic=[1:length(species)];
                jODE = [];
                timeFunName = 'timeDependentFactor';
                jntFactorName = 'jointDependentFactor';
            end

            n_reactions = size(stoichMatrix,2);
            n_pars = size(nonXTpars,1);
            speciesStoch = setdiff(species,upstreamODEs,'stable');

            varODEs = sym('varODEs',[1,length(upstreamODEs)]);

            % Delete previous propensity function m-files
            if ~exist([pwd,'/tmpPropensityFunctions'],'dir')
                mkdir([pwd,'/tmpPropensityFunctions'])
            end
            delete(append(pwd,'/tmpPropensityFunctions/',prefixName,'*'));
            addpath([pwd,'/tmpPropensityFunctions/'],'-begin')

            obj = cell(1,n_reactions);
            expr_t_vec = sym('w',[n_reactions,1]);
            expr_x_vec = sym('wx',[n_reactions,1]);
            expr_dt_vec_dode = sym('wx',[n_reactions,length(jODE)]);
            
            if computeSens
                expr_t_vec_sens = sym('w',[n_reactions,n_pars]);
                expr_x_vec_sens = sym('wx',[n_reactions,n_pars]);
            end

            anyLogical = zeros(1,n_reactions,'logical');
           
            t = sym('t','real');
            for iRxn = 1:n_reactions
                prop_vars = symvar(symbolicExpression{iRxn});
                syms(prop_vars,'real');
            end
            if ~isempty(upstreamODEs)
                syms(upstreamODEs,'positive')
            end
            oneSym = str2sym('1');

            % change to parfor?
            for iRxn = 1:n_reactions
                prop_vars = symvar(symbolicExpression{iRxn});
                hybridFactor =[];
                prefixNameLocal = [prefixName,'_',num2str(iRxn)];

                obj{iRxn} = ssit.Propensity(stoichMatrix(jStochastic,iRxn), iRxn);
                obj{iRxn}.ODEstoichVector = stoichMatrix(jODE,iRxn);

                if ~isempty(jODE)&&~isempty(obj{iRxn}.stoichVector)&&max(abs(obj{iRxn}.stoichVector))~=0&&max(abs(obj{iRxn}.ODEstoichVector))~=0
                    warning(['Reaction ',num2str(iRxn),' changes both ODE and stochastic species. Removing effect on upstream species.'])
                    obj{iRxn}.ODEstoichVector = 0*obj{iRxn}.ODEstoichVector;
                end
                % prop_vars = symvar(symbolicExpression{iRxn});

                % syms t real

                expr_tx = symbolicExpression{iRxn};
                % Convert to callable function
                Jx = zeros(1,length(prop_vars),'logical'); % List of variables corrresponding to species
                Jt = zeros(1,length(prop_vars),'logical'); % List of variables corrresponding to time
                if ~isempty(speciesStoch)
                    for j = 1:length(prop_vars)
                        Jx(j) = max(strcmp(string(prop_vars(j)),speciesStoch));
                        Jt(j) = max(strcmp(string(prop_vars(j)),['t',upstreamODEs]));
                    end
                end

                if sum(Jx)==0
                    % No x dependance
                    factors = [oneSym,symbolicExpression{iRxn}];
                elseif sum(Jt)==0
                    % No t dependence
                    factors = [symbolicExpression{iRxn},oneSym];
                else
                    % Determine if the symbolicExpression is isFactorizable
                    factors = factor(symbolicExpression{iRxn}, prop_vars(Jx));
                    factors = [prod(factors(2:end)),factors(1)];

                    % if ~isempty(upstreamODEs)
                    %     syms(upstreamODEs,'positive')
                    % end

                    if (ismember(t, symvar(factors(1))))||(~isempty(upstreamODEs)&&max(ismember(upstreamODEs,symvar(factors(1)))))
                        obj{iRxn}.isFactorizable = false;
                    end
                    for i2 = 2:length(factors)
                        fvars = symvar(factors(i2));
                        % TMP = string(fvars);
                        Jx = zeros(1,length(fvars),'logical');
                        for j = 1:length(fvars)
                            Jx(j) = max(strcmp(string(fvars(j)),speciesStoch));
                        end
                        % Jx = strcmp(TMP,speciesStoch);
                        if (sum(Jx) > 1)
                            obj{iRxn}.isFactorizable = false;
                            break;
                        end
                    end
                end

                expr_t = prod(factors(2:end));
                expr_x = factors(1);

                % Convert these symbolic expressions to anonymous
                % functions
                for i2=1:length(upstreamODEs)
                    expr_t = subs(expr_t,upstreamODEs{i2},varODEs(i2));
                    expr_x = subs(expr_x,upstreamODEs{i2},varODEs(i2));
                end

                if ~isempty(string(symvar(expr_t)))
                    if ~max(contains(string(symvar(expr_t)),'t'))&&~max(contains(string(symvar(expr_t)),'logT'))
                        % Check that there is actually a t-dependent reaction
                        % and otherwise combine.
                        expr_x = expr_x*expr_t;
                        expr_t=sym(1);
                    end
                end


                if (~isempty(string(symvar(expr_x)))&&max(contains(string(symvar(expr_x)),'logT')))||...
                        (~isempty(string(symvar(expr_t)))&&max(contains(string(symvar(expr_t)),'logX')))
                    obj{iRxn}.isFactorizable = false;
                    anyLogical(iRxn) = true;
                end

                if (obj{iRxn}.isFactorizable)
                    % If the propensity is isFactorizable, collect the
                    % time-dependent and state-dependent function handles


                    % The following commands will generate callable
                    % functions to compute propensity functions.  If there
                    % are no logical expressions, these are created as
                    % m-files, and if there are logical commands, these
                    % are created as much slower anonymous functions. this
                    % will require keeping track since the variables need
                    % to be sent individually to the anonymous functions.
                    if ~isempty(logicTerms{iRxn})&&(isfield(logicTerms{iRxn},'logT')||isfield(logicTerms{iRxn},'logJ'))
                        obj{iRxn}.anonymousT = true;
                        hybridFactor = sym2propfun(expr_t, true, false, nonXTpars(:,1), speciesStoch, varODEs, logicTerms(iRxn));
                        anyLogical(iRxn) = true;
                    % elseif sum(contains(string(symvar(expr_t)),'logX'))
                    %     obj{iRxn}.anonymousT = true;
                    %     obj{iRxn}.isFactorizable = false;
                    %     hybridFactor = sym2propfun(expr_t, true, false, nonXTpars(:,1), speciesStoch, varODEs, logicTerms(iRxn));
                    %     anyLogical(iRxn) = true;
                    else
                        obj{iRxn}.anonymousT = false;
                        % hybridFactor = sym2mFun(expr_t, true, false, nonXTpars(:,1), speciesStoch, varODEs);
                    end

                    if ~isempty(logicTerms{iRxn})&&(isfield(logicTerms{iRxn},'logX')||isfield(logicTerms{iRxn},'logJ'))
                        obj{iRxn}.anonymousX = true;
                        anyLogical(iRxn) = true;
                    % elseif sum(contains(string(symvar(expr_x)),'logT'))
                    %     obj{iRxn}.isFactorizable = false;
                    %     obj{iRxn}.anonymousX = true;
                    %     anyLogical(iRxn) = true;
                    else
                        obj{iRxn}.anonymousX = false;
                    end


                    % Compute the time-varying factor.
                    if ~isempty(logicTerms{iRxn})&&(isfield(logicTerms{iRxn},'logT')||isfield(logicTerms{iRxn},'logE'))
                        if ~isempty(jODE)
                            TmpHybridFactor = hybridFactor(rand,rand(size(nonXTpars(:,1))),rand(size(varODEs)));
                        else
                            TmpHybridFactor = hybridFactor(rand,rand(size(nonXTpars(:,1))));
                        end
                    else
                        if ~isempty(jODE)
                            obj{iRxn}.isTimeDependent = true;
                            % TmpHybridFactor =  hybridFactor(rand,[nonXTpars{:,2}]',rand(1,length(upstreamODEs)));
                            TmpHybridFactor = subs(expr_t,t,rand);
                            TmpHybridFactor = subs(TmpHybridFactor,nonXTpars(:,1),rand(size(nonXTpars(:,2))));
                            TmpHybridFactor = double(subs(TmpHybridFactor,varODEs,rand(size(varODEs))));
                        else
                            TmpHybridFactor = subs(expr_t,t,rand);
                            TmpHybridFactor = double(subs(TmpHybridFactor,nonXTpars(:,1),rand(size(nonXTpars(:,2)))));
                        end
                    end

                    % The time dependent signal must always be
                    % non-negative, so we check the sign and multiple by -1
                    % if necessary. 
                    % First for the time varying factors.
                    signHybridFactor = ((TmpHybridFactor>=0)-(TmpHybridFactor<0));
                    if obj{iRxn}.anonymousT
                        [obj{iRxn}.(timeFunName),expr_dt_vec_dodei] = ...
                            sym2propfun(signHybridFactor*expr_t, true, false, nonXTpars(:,1), speciesStoch, varODEs, logicTerms(iRxn), true);
                    else
                        [~,expr_dt_vec_dodei] = ...
                            sym2mFun(signHybridFactor*expr_t, true, false, nonXTpars(:,1), speciesStoch, varODEs, true, false, prefixNameLocal);
                    end

                    % Then for the state varying factors.
                    if obj{iRxn}.anonymousX
                        obj{iRxn}.stateDependentFactor =...
                            sym2propfun(signHybridFactor*expr_x, false, true, nonXTpars(:,1), speciesStoch, [], logicTerms(iRxn));                    
                    else
                        obj{iRxn}.stateDependentFactor =...
                            sym2mFun(signHybridFactor*expr_x, false, true, nonXTpars(:,1), speciesStoch, [], false, true, prefixNameLocal);
                    end

                    expr_t_vec(iRxn) = signHybridFactor*expr_t;
                    expr_x_vec(iRxn) = signHybridFactor*expr_x;
                    if ~isempty(expr_dt_vec_dodei)
                        expr_dt_vec_dode(iRxn,:) = signHybridFactor*expr_dt_vec_dodei;
                    end
                    if ismember(t,symvar(expr_t))||isfield(logicTerms{iRxn},'logT')&&~isempty(logicTerms{iRxn}.logT)
                        obj{iRxn}.isTimeDependent = true;
                    end
                    if computeSens
                        for iPar = 1:n_pars
                            expr_t_vec_sens(iRxn,iPar) = signHybridFactor*diff(expr_t,nonXTpars{iPar,1});
                            expr_x_vec_sens(iRxn,iPar) = signHybridFactor*diff(expr_x,nonXTpars{iPar,1});
                        end
                    end

                else
                    for i2=1:length(upstreamODEs)
                        expr_tx = subs(expr_tx,upstreamODEs{i2},varODEs(i2));
                    end

                    [obj{iRxn}.(jntFactorName),expr_dt_vec_dodei] = ...
                        sym2propfun(expr_tx, true, true, nonXTpars(:,1), speciesStoch, varODEs, logicTerms(iRxn), true);
                    obj{iRxn}.isTimeDependent = true;
                    expr_t_vec(iRxn) = sym(NaN);
                    expr_x_vec(iRxn) = sym(NaN);
                    
                    if ~isempty(expr_dt_vec_dodei)
                        expr_dt_vec_dode(iRxn,:) = expr_dt_vec_dodei;
                    end
                    
                    if computeSens
                        for iPar = 1:n_pars
                            expr_tx_vec_sens(iRxn,iPar) = diff(expr_tx,nonXTpars{iPar,1});
                        end                   
                    end

                end
            end
            
            prefixNameLocal = prefixName;

            % Create a vector function for the time varying part of the
            % propensity functions so that these can be found all at once.
            if any(anyLogical)
                hybridFactorVector = sym2propfun(expr_t_vec, true, false, nonXTpars(:,1), speciesStoch, varODEs, logicTerms);
                xFactorVector = sym2propfun(expr_x_vec, false, true, nonXTpars(:,1), speciesStoch, varODEs, logicTerms);
                if ~isempty(expr_dt_vec_dode)
                    obj{1}.DhybridFactorDodesVec = sym2propfun(expr_dt_vec_dode, true, false, nonXTpars(:,1), speciesStoch, varODEs, logicTerms);
                end
                if computeSens
                    obj{1}.sensStateFactorVec = sym2propfun(expr_x_vec_sens, false, true, nonXTpars(:,1), speciesStoch, varODEs, logicTerms);
                    for iRxn = 1:n_reactions
                        obj{iRxn}.sensStateFactor = cell(1,n_pars);
                        for ipar = 1:n_pars
                            obj{iRxn}.sensStateFactor{ipar} =  sym2propfun(expr_x_vec_sens(iRxn,ipar), false, true, nonXTpars(:,1), speciesStoch, varODEs, logicTerms);
                            obj{1}.sensTimeFactorVec{ipar} = sym2propfun(expr_t_vec_sens(:,ipar), true, false, nonXTpars(:,1), speciesStoch, varODEs, logicTerms);
                        end
                    end
                end
            else
                hybridFactorVector = sym2mFun(expr_t_vec, true, false, nonXTpars(:,1), speciesStoch, varODEs, false, true, [prefixNameLocal,'_t']);
                xFactorVector = sym2mFun(expr_x_vec, false, true, nonXTpars(:,1), speciesStoch, varODEs, false, true, [prefixNameLocal,'_x']);
                if ~isempty(expr_dt_vec_dode)
                    obj{1}.DhybridFactorDodesVec = sym2mFun(expr_dt_vec_dode, true, false, nonXTpars(:,1), speciesStoch, varODEs, false, true, [prefixNameLocal,'_dt']);
                end
                if computeSens
                    obj{1}.sensStateFactorVec = sym2mFun(expr_x_vec_sens, false, true, nonXTpars(:,1), speciesStoch, varODEs, false, true, [prefixNameLocal,'_s']);
                    poolobj = gcp("nocreate");
                    if ~isempty(poolobj)&&n_reactions>0
                        parfor iRxn = 1:n_reactions
                            obj{iRxn}.sensStateFactor = cell(1,n_pars);
                            for ipar = 1:n_pars
                                prefixNameLocal = [prefixName,'_',num2str(iRxn),'_',num2str(ipar)];
                                obj{iRxn}.sensStateFactor{ipar} =  sym2mFun(expr_x_vec_sens(iRxn,ipar), false, true, nonXTpars(:,1), speciesStoch, varODEs, false, true, prefixNameLocal);
                            end
                        end
                    else
                        for iRxn = 1:n_reactions
                            obj{iRxn}.sensStateFactor = cell(1,n_pars);
                            for ipar = 1:n_pars
                                prefixNameLocal = [prefixName,'_',num2str(iRxn),'_',num2str(ipar)];
                                obj{iRxn}.sensStateFactor{ipar} =  sym2mFun(expr_x_vec_sens(iRxn,ipar), false, true, nonXTpars(:,1), speciesStoch, varODEs, false, true, prefixNameLocal);
                            end
                        end
                    end
                    for ipar = 1:n_pars
                        prefixNameLocal = [prefixName,'_v_',num2str(1),'_',num2str(ipar)];
                        obj{1}.sensTimeFactorVec{ipar} = sym2mFun(expr_t_vec_sens(:,ipar), true, false, nonXTpars(:,1), speciesStoch, varODEs, false, true, prefixNameLocal);
                    end
                end
            end
            obj{1}.hybridFactorVector = hybridFactorVector;
            obj{1}.xFactorVector = xFactorVector;
        end

        function obj = createFromString(fstr, stoichVector, reactionIndex)
            % Construct an instance of Propensity from a symbolic expression.
            %
            % Parameters
            % ----------
            %
            %   fstr: string
            %       string expression of the propensity.
            %
            %   stoichVector: column vector
            %       is the column stoichiometry vector associated
            %       with this propensity.
            %
            %   reactionIndex: integer
            %       index of the reaction associated with this propensity.
            %
            % Returns
            % -------
            %
            %   obj: Propensity
            %       an instance of the Propensity class.
            %
            import ssit.*
            fstr = strrep(fstr,'..','.');
            [ft, fx, isFactorizable] = separateExpression(fstr);

            if isFactorizable
                ft = strrep(ft,'..','.');fx = strrep(fx,'..','.');
                obj = Propensity.createAsSeparable(str2func(['@(t)' ft]), str2func(['@(x)' fx]),...
                    stoichVector, reactionIndex);
            else
                obj = Propensity.createAsNotSeparable(str2func(['@(t,x)' fstr]), stoichVector, reactionIndex, fstr);
            end
        end

        function obj = createAsSeparable(fthandle, fxhandle, stoichVector, reactionIndex)
            % Create an object for a time-varying propensity that is
            % separable.
            %
            % Parameters
            % ----------
            %
            %   fthandle: function handle
            %       time-depdendent factor of the propensity function.
            %
            %   fxhandle: function handle
            %       state-dependent factor of the propensity function.
            %
            %   stoichVector: column vector
            %       stoichiometry vector associated with this propensity.
            %
            %   reactionIndex: integer
            %       index of the reaction associated with this propensity.
            %
            % Returns
            % -------
            %
            %   obj: Propensity
            %       an instance of Propensity, where obj.isFactorizable ==
            %       true.
            %
            import ssit.*
            obj = Propensity(stoichVector, reactionIndex);
            if strcmp(func2str(fthandle),'@(t)1')
                obj.isTimeDependent = false;
            else
                obj.isTimeDependent = true;
            end
            obj.isFactorizable = true;
            obj.timeDependentFactor = fthandle;
            obj.stateDependentFactor = fxhandle;
        end
    end

    methods (Static)
        function [stNew,logicTerms,counter] = stripLogicals(st,species,counter)
            logTypes = {'=','>','<'};
            logicTerms = [];
            n = [0,0,0];
            stNew = st;
            for i=1:3
                while ~isempty(strfind(stNew,logTypes{i}))
                    J = strfind(stNew,logTypes{i});
                    for j = 1%:length(J)
                        logvals = (stNew == '(') - (stNew == ')');
                        k1 = J(j);
                        pos = 0;
                        while pos<1
                            k1 = k1-1;
                            pos = sum(logvals(k1:J(j)));
                        end
                        k2 = J(j);
                        neg = 0;
                        while neg>-1
                            k2 = k2+1;
                            neg = sum(logvals(J(j):k2));
                        end

                        % K = strfind(stNew,'(');
                        % k1 = max(K(K<J(j)));
                        % K = strfind(stNew,')');
                        % k2 = min(K(K>J(j)));
                        logE = stNew(k1:k2);
                        if ~isempty(regexp(logE,'\<t\>'))&&max(contains(logE,species))
                            n(1)=n(1)+1;
                            logicTerms.logJ{n(1),1} = logE;
                            counter = counter+1;
                            logicTerms.logJ{n(1),2} = ['logJ',num2str(counter)];
                            stNew = strrep(stNew,logE,['(',logicTerms.logJ{n(1),2},')']);
                            % stNew = regexprep(stNew,['\<',logE,'\>'],['(',logicTerms.logJ{n(1),2},')']);
                        elseif ~isempty(regexp(logE,'\<t\>'))
                            n(2)=n(2)+1;
                            logicTerms.logT{n(2),1} = logE;
                            counter = counter+1;
                            logicTerms.logT{n(2),2} = ['logT',num2str(counter)];
                            stNew = strrep(stNew,logE,['(',logicTerms.logT{n(2),2},')']);
                            % stNew = regexprep(stNew,['\<',logE,'\>'],['(',logicTerms.logT{n(2),2},')']);
                        elseif max(contains(logE,species))
                            n(3)=n(3)+1;
                            logicTerms.logX{n(3),1} = logE;
                            counter = counter+1;
                            logicTerms.logX{n(3),2} = ['logX',num2str(counter)];
                            stNew = strrep(stNew,logE,['(',logicTerms.logX{n(3),2},')']);
                            % stNew = regexprep(stNew,['\<',logE,'\>'],['(',logicTerms.logX{n(3),2},')']);
                        end
                    end
                end
            end
        end
    end
end

function [ft, fx, isFactorizable] = separateExpression(expr)
% cft = [];
fx = [];
isFactorizable = [];
% Check input
left_bracket_loc = strfind(expr, '{');
right_bracket_loc = strfind(expr, '}');

if (length(left_bracket_loc) > length(right_bracket_loc))
    fprintf('Syntax error, missing ''}'' in one of the propensities.\n');
    return;
elseif (length(left_bracket_loc) < length(right_bracket_loc))
    fprintf('Syntax error, missing ''{'' in one of the propensities.\n');
    return;
end

if (length(left_bracket_loc) > 2)
    fprintf('We currently only support separating propensity into two factor.\n');
    return;
end

if (isempty(left_bracket_loc))
    % Is it a function of only t?
    if (~contains(expr, 't'))
        isFactorizable = true;
        ft = '1';
        fx = expr;
        return;
    elseif (~contains(regexprep(expr,'\<exp\>','EP'), 'x')) % Check if it has 'x' and not 'exp'
        isFactorizable = true;
        ft = expr;
        fx = '1';
        return;
    else
        isFactorizable = false;
        ft = [];
        fx = [];
        return;
    end
end
%
isFactorizable = true;
f1 = expr(left_bracket_loc(1)+1:right_bracket_loc(1)-1);
f2 = expr(left_bracket_loc(2)+1:right_bracket_loc(2)-1);

if (contains(f1, 't') || contains(f1, 'I'))
    ft = f1;
    fx = f2;
else
    ft = f2;
    fx = f1;
end
end

function [exprHandle,exprJac] = sym2propfun(symbolicExpression, time_dep, state_dep, nonXTpars, species, varODEs, logicTerms, jacobian)
arguments
    symbolicExpression
    time_dep
    state_dep
    nonXTpars
    species
    varODEs = []
    logicTerms = {}
    jacobian = false
end
% Convert a symbolic expression into usable anonoymous function
% handle.
% We assume that the time variable is always named 't' and the
% state variables are always named with the form 'x*' where the
% expression following the character 'x' tells us the index of
% the species. For example, 'x1', 'x2', 'x12' etc.

% To do:  Adjust this to first factorize the symbolic expression. Check
% to see if terms with x1, x2, etc can be separated completely from
% those with t or other veriables.  If so, lump all these non-xi terms together
% in an anonymous functon (of all non-x parameters) and put all the
% x-dependent terms into another function. Will need to check that this
% does not disrupt codes for sensitivity analysis.

% import ssit.fsp.*
varNames = string(symvar(symbolicExpression));
varNames = unique([varNames,species{:}]);
% len = zeros(1,length(varNames));
% for i = 1:length(varNames)
%     len(i) = length(varNames{i});
% end
% [~,J] = sort(len,'descend');
% varNames = varNames(J);
% sort(varNames, 'descend');

if jacobian&&~isempty(varODEs)
    exprJac = sym('dedx',[1,length(varODEs)]);
    for i=1:length(varODEs)
        exprJac(i) = diff(symbolicExpression,varODEs(i));
    end
else
    exprJac=[];
end


exprStr = char(symbolicExpression);

% Get rid of  max rules.
% k = strfind(exprStr,', ''omitnan');
% if ~isempty(k)
%     exprStr = [exprStr(1:k-1),')'];
% end
exprStr = strrep(exprStr,', ''omitnan'', false','');
exprStr = strrep(exprStr,', \"''omitnan''\", false','');
% Same but for another format.
% k = strfind(exprStr,', \"''omitnan''\"');
% if ~isempty(k)
%     exprStr = [exprStr(1:k-1),')'];
% end

% exprStr = char(strrep(exprStr,", [], 2, 'omitnan', false",""));
opVar = {'*','/','^'};
for i = 1:length(opVar)
    op = opVar{i};
    exprStr = strrep(exprStr, op, ['.' op]);
end

% parStr = [];
% for iPar = 1:length(nonXTpars)
%     parStr = [parStr,', ',nonXTpars{iPar}];
% end
parStr = ', Parameters';

if ~isempty(varODEs)
    parStr = [parStr,', varODEs'];
    for i=length(varODEs):-1:1
        % exprStr = strrep(exprStr, ['varODEs',num2str(i)], ['varODEs(',num2str(i),')']);
        exprStr = regexprep(exprStr, ['\<varODEs',num2str(i),'\>'], ['varODEs(',num2str(i),')']);
    end
end

for i=1:length(logicTerms)
    if isfield(logicTerms{i},'logT')
        for j=1:size(logicTerms{i}.logT,1)
            % exprStr=strrep(exprStr,logicTerms{i}.logT{j,2},logicTerms{i}.logT{j,1});
            exprStr=regexprep(exprStr,['\<',logicTerms{i}.logT{j,2},'\>'],logicTerms{i}.logT{j,1});
        end
    end
    if isfield(logicTerms{i},'logX')
        %             state_dep = true;
        for j=1:size(logicTerms{i}.logX,1)
            % exprStr=strrep(exprStr,logicTerms{i}.logX{j,2},logicTerms{i}.logX{j,1});
            exprStr=regexprep(exprStr,['\<',logicTerms{i}.logX{j,2},'\>'],logicTerms{i}.logX{j,1});
        end
    end
    if isfield(logicTerms{i},'logJ')
        %             time_dep = true;
        %             state_dep = true;
        for j=1:size(logicTerms{i}.logX,1)
            % exprStr=strrep(exprStr,logicTerms{i}.logX{j,2},logicTerms{i}.logX{j,1});
            exprStr=regexprep(exprStr,['\<',logicTerms{i}.logX{j,2},'\>'],logicTerms{i}.logX{j,1});
        end
    end
end

lens = zeros(1,length(nonXTpars));
for i = 1:length(nonXTpars)
    lens(i) = length(nonXTpars{i});
end
[~,J] = sort(lens,'descend');

for i = 1:length(nonXTpars)
    % exprStr = strrep(exprStr, nonXTpars{J(i)}, ['$(',num2str(J(i)),')']);
    exprStr = regexprep(exprStr, ['\<',nonXTpars{J(i)},'\>'], ['Parameters(',num2str(J(i)),')']);
end
% exprStr = strrep(exprStr,'$','Parameters');

for i = length(nonXTpars):-1:1
    exprStr = regexprep(exprStr, ['\<parameters',num2str(i),'\>'], ['Parameters(',num2str(i),')']);
end


if (time_dep && state_dep)
    fhandle_var = ['@(t, x',parStr,')'];
    for i = 1:length(varNames)
        old_name = char(varNames(i));
        if max(ismember(species, old_name))
            j = find(ismember(species, old_name));
            new_name = ['x(', num2str(j), ', :)'];
            % exprStr = strrep(exprStr, old_name, new_name);
            exprStr = regexprep(exprStr, ['\<',old_name,'\>'], new_name);
        end
    end
    exprHandle = str2func([fhandle_var exprStr]);

elseif (time_dep)
    fhandle_var = ['@(t',parStr,')'];
    exprHandle = str2func([fhandle_var '(' exprStr ')']);
else
    fhandle_var = ['@(x',parStr,')'];
    if ~isempty(varNames)
        for i = 1:length(varNames)
            old_name = char(varNames(i));
            if max(ismember(species, old_name))
                j = find(ismember(species, old_name));
                new_name = ['x(', num2str(j), ', :)'];
                % exprStr = strrep(exprStr, old_name, new_name);
                exprStr = regexprep(exprStr, ['\<',old_name,'\>'], new_name);
            end
        end
    else
        exprStr = ['( ' exprStr ').* ones(1, size(x,2))'];
    end
    exprHandle = str2func([fhandle_var exprStr]);
end
end
function [exprHandle,exprJac] = sym2mFun(symbolicExpression, time_dep, state_dep, nonXTpars, species, varODEs, jacobian, writeFiles, prefixName)
arguments
    symbolicExpression
    time_dep
    state_dep
    nonXTpars
    species
    varODEs = {}
    jacobian = false
    writeFiles = true
    prefixName = []
end
% This function writes an executable m-file for the provided expression.
varNames = string(symvar(symbolicExpression));
states = sym("states",[length(species),1],'positive');
parameters = sym("parameters",[length(nonXTpars),1],'positive');
syms t real
for i = 1:length(varNames)
    old_name = char(varNames(i));
    if max(ismember(species, old_name))
        j = find(ismember(species, old_name));
        symbolicExpression = subs(symbolicExpression,species{j},states(j));
    end
    if max(ismember(nonXTpars, old_name))
        j = find(ismember(nonXTpars, old_name));
        symbolicExpression = subs(symbolicExpression,nonXTpars{j},parameters(j));
    end
end

if isempty(prefixName)
    prefixName = pwd;
    j = find(prefixName==filesep,1,'last');
    prefixName = prefixName(j+1:end);
end

ifn = sum(contains({dir('tmpPropensityFunctions').name},[prefixName,'_fun']))+1;
fn = [pwd,'/tmpPropensityFunctions/',prefixName,'_fun_',num2str(ifn),'.m'];

if jacobian&&~isempty(varODEs)
    exprJac = sym(zeros([1,length(varODEs)]));
    [~,a2] = intersect(varODEs,symvar(symbolicExpression));
    for i=1:length(a2)
        exprJac(a2(i)) = diff(symbolicExpression,varODEs(a2(i)));
    end
else
    exprJac=[];
end

if writeFiles
    if isempty(varODEs)
        if (time_dep && state_dep)
            exprHandle = matlabFunction(symbolicExpression,'Vars',{t,states,parameters},'File',fn,'Sparse',true);
        elseif time_dep
            exprHandle = matlabFunction(symbolicExpression,'Vars',{t,parameters},'File',fn,'Sparse',true);
        else
            exprHandle = matlabFunction(symbolicExpression,'Vars',{states,parameters},'File',fn,'Sparse',false);
        end
    else
        if (time_dep && state_dep)
            exprHandle = matlabFunction(symbolicExpression,'Vars',{t,states,parameters,varODEs},'File',fn,'Sparse',true);
        elseif time_dep
            exprHandle = matlabFunction(symbolicExpression,'Vars',{t,parameters,varODEs},'File',fn,'Sparse',true);
        else
            exprHandle = matlabFunction(symbolicExpression,'Vars',{states,parameters,varODEs},'File',fn,'Sparse',true);
        end
    end
    if ~strcmp(which(func2str(exprHandle)),fn)
        disp('WARNING -- it appears that new propensity functions is redundant to one already on the search path.')
        disp(['at: ',which(func2str(exprHandle))]);
    end
else
    exprHandle=[];
end
% if jacobian
%     fn = [pwd,'/tmpPropensityFunctions/',prefix,'_fun_',num2str(ifn+1),'.m'];
%     if isempty(varODEs)
%         % if (time_dep && state_dep)
%         %     exprHandleJac = matlabFunction(exprJac,'Vars',{t,states,parameters},'File',fn);
%         % elseif state_dep
%         %     exprHandleJac = matlabFunction(exprJac,'Vars',{states,parameters},'File',fn);
%         % end
%     else
%         if (time_dep && state_dep)
%             exprHandleJac = matlabFunction(exprJac,'Vars',{t,states,parameters,varODEs},'File',fn);
%         elseif time_dep
%             exprHandleJac = matlabFunction(exprJac,'Vars',{states,parameters,varODEs},'File',fn);
%         end
%     end
% end
end
