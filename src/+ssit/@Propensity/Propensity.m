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
        originalString;
    end

    methods (Access = public)
        function obj = Propensity(stoichVector, reactionIndex)
            obj.reactionIndex = reactionIndex;
            obj.stoichVector = stoichVector;
            obj.isTimeDependent = false;
            obj.isFactorizable = true;
        end

        function y = evaluate(obj, t, x)
            if (size(x, 1) ~= size(obj.stoichVector, 1))
                error('Propensity.eval : input state dimension is not consistent with the stoichiometry vector in the object''s record');
            end
            if (~obj.isTimeDependent)
                y = obj.stateDependentFactor(x);
            else
                if (obj.isFactorizable)
                    y = obj.timeDependentFactor(t)*obj.stateDependentFactor(x);
                else
                    try
                        y = obj.jointDependentFactor(t, x);
                    catch
%                         f = uifigure;
%                         X = uiconfirm(f, 'Differentiation is not possible for your current propensity functions. Please change to finite difference sensitivity, or try with pieceswise functions.', 'Non-Differentiable',...
%                         'Options', {'Okay - Quit','Try with Piecewise Functions'},...
%                         'DefaultOption','Okay - Quit');
%                         close(f)
%                         if strcmp(X,'Okay - Quit')
%                             error('Sensitivity not analytically computable in current version. Please use Finite Difference');
%                         end
%                         disp('trying with piecewise funtions');
                        TMP = func2str(obj.jointDependentFactor);
                        TMP = strrep(TMP,'t','tt');
                        TMP = strrep(TMP,'x(1,:)','x1');
                        TMP = strrep(TMP,'x(2,:)','x2');
                        TMP = strrep(TMP,'x(3,:)','x3');
                        syms tt
                        syms x1 x2 x3 positive
                        eval(['Q=',TMP(8:end),';'])
                        y = eval(subs(Q,{tt,x1,x2,x3},{t*ones(1,size(x,2)),x(1,:),x(2,:),x(3,:)}));
                    end
                    
                end
            end
        end

        function y = evaluateStateFactor(obj, x)
            if (size(x, 1) ~= size(obj.stoichVector, 1))
                error('Propensity.evalAt : input state dimension is not consistent with the stoichiometry vector in the object''s record');
            end
             try
                y = obj.stateDependentFactor(x);
            catch
                TMP = func2str(obj.stateDependentFactor);
                TMP = strrep(TMP,'t','tt');
                TMP = strrep(TMP,'x(1,:)','x1');
                TMP = strrep(TMP,'x(2,:)','x2');
                TMP = strrep(TMP,'x(3,:)','x3');
                TMP = strrep(TMP,'x(4,:)','x4');
                syms tt
                syms x1 x2 x3 x4 positive
                eval(['Q=',TMP(5:end),';'])
                y = eval(subs(Q,{x1,x2,x3,x4},{x(1,:),x(2,:),x(3,:),x(4,:)}));
            end
  
        end
    end
    methods (Access = public, Static)
        function obj = createFromSym(symbolicExpression, stoichVector, reactionIndex)
            % Construct an instance of Propensity from a symbolic expression.
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
            % Returns
            % -------
            %
            %   obj: Propensity
            %       an instance of the Propensity class.
            %
            obj = ssit.Propensity(stoichVector, reactionIndex);


            % Now determine time-dependency
            prop_vars = symvar(symbolicExpression);
            str_tmp = string(prop_vars);
            for i = 1:length(str_tmp)
                if (strcmp(str_tmp(i), "t"))&&~(strcmp(str_tmp, "@(t)1"))
                    obj.isTimeDependent = true;
                    break;
                end
            end

            if (obj.isTimeDependent)
                % Determine if the symbolicExpression is isFactorizable
                syms t;
                factors = factor(symbolicExpression, t);
                if (length(factors) == 1)
                    factors = [str2sym('1') factors];
                end
                if (ismember(t, symvar(factors(1))))
                    obj.isFactorizable = false;
                end
                for i = 2:length(factors)
                    if (length(symvar(factors(i))) > 1)
                        obj.isFactorizable = false;
                        break;
                    end
                end

                if (obj.isFactorizable)
                    % If the propensity is isFactorizable, collect the
                    % time-dependent and state-dependent function handles
                    expr_t = prod(factors(2:end));
                    expr_x = factors(1);
                    % Convert these symbolic expressions to anonymous
                    % functions
                    obj.timeDependentFactor = sym2propfun(expr_t, true, false);
                    obj.stateDependentFactor = sym2propfun(expr_x, false, true);
                else
                    % If it is not factorizable, collect everything into a
                    % single function handle
                    expr_tx = symbolicExpression;
                    % Convert to anonymous function
                    obj.jointDependentFactor = sym2propfun(expr_tx, true, true);
                end
            else
                expr_x = symbolicExpression;
                % Convert to anonymous function
                obj.timeDependentFactor = @(t) 1;
                obj.stateDependentFactor = sym2propfun(expr_x, false, true);
            end
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

        function obj = createAsTimeInvariant(fxhandle, stoichVector, reactionIndex)
            % Create an object for a time-invariant propensity function.
            %
            % Parameters
            % ----------
            %
            %   fxhandle: function handle
            %       state-dependent factor of the propensity function. Must have the same syntax as the property Propensity.stateDependentFactor.
            %
            %   stoichVector: column vector
            %       stoichiometry vector associated with this propensity.
            %
            %   reactionIndex: integer
            %       index of the reaction associated with this propensity.

            import ssit.*
            obj = Propensity(stoichVector, reactionIndex);
            obj.isTimeDependent = false;
            obj.isFactorizable = true;
            obj.stateDependentFactor = fxhandle;
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

        function obj = createAsNotSeparable(ftxhandle, stoichVector, reactionIndex, str)
            % Create an object for a time-varying propensity that is cannot
            % be factorized into a time-varying factor and a time-invariant
            % factor.
            %
            % Parameters
            % ----------
            %
            %   ftxhandle: function handle
            %       this function handle must evaluate the propensity at
            %       a given joint value of the time and state variables.
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
            %       false.
            %
            import ssit.*
            obj = Propensity(stoichVector, reactionIndex);
            obj.isTimeDependent = true;
            obj.isFactorizable = false;
            obj.jointDependentFactor = ftxhandle;
        end
    end
end

function [ft, fx, isFactorizable] = separateExpression(expr)
    cft = [];
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
        % Is it a function of only x?
        if (~contains(expr, 't'))
            isFactorizable = true;
            ft = '1';
            fx = expr;
            return;
        elseif (~contains(expr, 'x'))
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
    %     ft = strip(ft, 'left', '{');
    %     ft = strip(ft, 'right', '}');
    %     fx = strip(fx, 'left', '{');
    %     fx = strip(fx, 'right', '}');
end

function y = sym2propfun(symbolicExpression, time_dep, state_dep)
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

    import ssit.fsp.*
    exprStr = char(symbolicExpression);
    varNames = string(symvar(symbolicExpression));
    sort(varNames, 'descend');
    opVar = {'*','/','^'};
    for i = 1:length(opVar)
        op = opVar{i};
        exprStr = strrep(exprStr, op, ['.' op]);
    end

    if (time_dep && state_dep)
        fhandle_var = '@(t, x)';
        for i = 1:length(varNames)
            old_name = char(varNames(i));
            if (ismember('x', old_name))
                new_name = [old_name(1) '(' old_name(2:end) ', :)'];
                exprStr = strrep(exprStr, old_name, new_name);
            end
        end
        y = str2func([fhandle_var exprStr]);
    elseif (time_dep)
        fhandle_var = '@(t) ';
        y = str2func([fhandle_var '(' exprStr ')']);
    else
        fhandle_var = '@(x)';
        if ~isempty(varNames)
            for i = 1:length(varNames)
                old_name = char(varNames(i));
                if (ismember('x', old_name))
                    new_name = [old_name(1) '(' old_name(2:end) ', :)'];
                    exprStr = strrep(exprStr, old_name, new_name);
                end
            end
        else
            exprStr = ['( ' exprStr ').* ones(1, size(x,2))'];
        end
        y = str2func([fhandle_var exprStr]);
    end
end
