classdef (HandleCompatible) SrnModel
    %SRNMODEL Storing information for stochastic reaction network
    %models whose propensities depend on parameters.
    %
    % Parameters
    % ----------
    %
    %   stoichiometry: NxM matrix
    %       Matrix of stoichiometry vectors. N is the number of species; M is the number of reactions.
    %
    %   parameterNames: cell array of strings
    %       Names of the model parameters.
    %
    %   propensityExpressions: cell array of strings
    %       Mathematical expressions for reaction propensities. This cell
    %       must have the same length as the number of reactions in the
    %       model. These strings can contain variables included in
    %       ```parameterNames``` and time-varying functions defined in
    %       ```timeVaryingInputExpressions```.
    %
    %   timeVaryingInputExpressions: 1-D cell array of strings
    %       Mathematical expressions for the time-varying functions.
    %
    %   propensityDerivativeExpressions: 2-D cell array of strings, optional
    %       Mathematical expressions for the partial derivatives of the
    %       propnesities with respect to the parameters. Number of columns
    %       must equal the number of parameters, and number of rows must
    %       equal the number of propensity functions. Users can leave the cell
    %       element empty if the derivative is zero.
    %
    properties
        stoichiometry
        parameterNames
        propensityExpressions
        timeVaryingInputExpressions
        propensityDerivativeExpressions
    end

    methods
        function obj = SrnModel(stoichiometry,...
                propensityExpressions,...
                parameterNames, ...
                inputExpressions, ...
                propensityDerivativeExpressions)
            %Construct an `SrnModel` instance.
            %   
            % Parameters
            % ----------
            %   
            %   stoichiometry: matrix
            %       Net stoichiometry matrix.
            %
            %   propensityExpressions: cell of strings
            %       Expressions for the propensity functions.
            %
            %   parameterNames: cell of strings
            %       Name of model parameters.
            %
            %   inputExpressions: cell of strings
            %       Expressions for the time-varying inputs used in
            %       `propensityExpressions`.
            %
            %   propensityDerivativeExpressions: 2-D cell array of strings
            %       

            arguments
                stoichiometry (:,:) {mustBeInteger}
                propensityExpressions cell
                parameterNames cell
                inputExpressions (:,:) cell=cell(0,0)
                propensityDerivativeExpressions (:,:) cell=cell(0,0)
            end

            obj.stoichiometry = stoichiometry;
            obj.propensityExpressions = propensityExpressions;
            obj.parameterNames = parameterNames;
            obj.timeVaryingInputExpressions = inputExpressions;
            obj.propensityDerivativeExpressions = propensityDerivativeExpressions;
        end

        %------------------------------------------------------------------
        function props = createPropensities(obj, parsDict, varNames)
            %CREATEPROPENSITY Create an instance of Propensity class.
            %   prop = obj.createPropensities(parameters)

            arguments
                obj
                parsDict
                varNames = {'x1';'x2';'x3'}
            end

            if isa(parsDict, 'cell')
                parsDict = containers.Map(parsDict(:,1), parsDict(:,2));
            end

            numProps = length(obj.propensityExpressions);
            propstrings = obj.processPropensityStrings(obj.propensityExpressions, ...
                obj.timeVaryingInputExpressions,...
                parsDict,...
                'str',varNames);

            props = cell(numProps, 1);
            for i= 1:numProps
                propstrings{i} = strrep(propstrings{i},'..','.');
                props{i} = ssit.Propensity.createFromString(propstrings{i},...
                    obj.stoichiometry(:,i),...
                    i);
            end
        end

        %------------------------------------------------------------------
        function [propensityDerivatives,isSensComputable] =...
                findPropensityDerivativesSymbolic(obj, parsDict, varNames)
            arguments
                obj
                parsDict
                varNames = {'x1';'x2';'x3'}
            end
 
            propensityCount = length(obj.propensityExpressions);
            parameterCount = length(obj.parameterNames);

            % Process the propensity strings to find their partial derivatives wrt parameters
            dwstr = cell(propensityCount, parameterCount);

            for i = 1:propensityCount
                tmp = obj.propensityExpressions{i};
                % Substitute the inputs for their expressions
                for j = 1:size(obj.timeVaryingInputExpressions,1)
                    nm = obj.timeVaryingInputExpressions{j,1};
                    va = obj.timeVaryingInputExpressions{j,2};
                    tmp = strrep(tmp, nm, ['(' va ')']);
                end
                
                 [dfk, ftx, isFactorizable] = separateExpressionParVxt(tmp,obj.parameterNames,varNames);
                 for j=1:length(isFactorizable)
                    if isFactorizable(j)
                        dwstr(i, j) = {(dfk{j}*ftx{j})};
                        if ~isempty(dwstr(i, j))
                            dwstr(i, j) = {char(dwstr{i, j})};
                        end
%                         dwstr2(i, j) = findDerivativesSyms(tmp, obj.parameterNames(j));
                    else
                        dwstr(i, j) = findDerivativesSyms(tmp, obj.parameterNames(j),varNames);
                    end
                 end
            end
                
            isSensComputable = ones(1,parameterCount);
            % Put parameter values to all function strings           
            for i = 1:propensityCount
                for j = 1:parameterCount
                    try
                        if ~isempty(dwstr{i, j})
                            dwstr(i, j) = obj.processPropensityStrings(dwstr(i, j), [], parsDict, "str", varNames);
                        end
                    catch
                        disp(['Propensity ',obj.parameterNames{j},' is not symbolically differentiable. Please try finite difference approach.'])
                    end
                end
            end

            propensityDerivatives = cell(propensityCount, parameterCount);
            for i = 1:propensityCount
                for j = 1:parameterCount
%                     if isSensComputable(j)
                        try
                            if ~isempty(dwstr{i,j})
                                propensityDerivatives{i,j} = ...
                                    ssit.Propensity.createFromString(dwstr{i,j}, obj.stoichiometry(:,i), i);
                            end
                        catch
                        disp(['Propensity ',obj.parameterNames{j},' is not computable. May give an error. Please try finite difference approach.'])
                        isSensComputable(j) = 0;
                        end
%                     end
                end
            end
        end
        %------------------------------------------------------------------
    end

    methods (Static)
        function props = processPropensityStrings(propensityExpressions,inputs,parDict,str_or_fun,varNames)
            % Convert the propensity string expressions into vectorized matlab functions. 
            %
            % Parameters
            % ----------
            %
            %   propensityExpressions: cell array of char arrays
            %       Here, propensities{j} is the expression for the propensity function
            %       of the j-th reaction.
            %
            %   inputs: cell array of char arrays
            %       Cell array of time-varying input formulas.
            %
            %   parDict: containers.Map object
            %       The dictionary of parameters, where each key is the name of a
            %       parameter (as a char array), and the associated value is the numerical value of that
            %       parameter.
            %
            %   str_or_fun: string
            %       If 'string', returns a char array expression for the formula. If 'fun',
            %       returns an anonoymous function.
            %
            % Returns
            % -------
            %
            %   props: cell array
            %       Each cell is the processed string expression for each propensity
            %       function (if str_or_fun is 'str' or unspecified), or an anoynoumous
            %       function (if str_or_fun is 'fun').
            %
            %
            % REMARKS: This function is confusing and it is serving two
            % purpose. The 'str' option gives us the vectorized
            % propensities for FSP, while the 'func' option gives us the
            % scalar propensities for SSA. It is better to break this into
            % two functions, each serve a single purpose (either for FSP or
            % for SSA).

            arguments
                propensityExpressions (:,1) cell
                inputs (:,2) cell={}
                parDict=[]
                str_or_fun char="str"
                varNames = {'x1','x2','x3'}
            end

            if ~isempty(parDict) && ~isa(parDict, 'containers.Map')
                parDict = containers.Map(parDict(:,1), parDict(:,2));
            end

            props = cell(length(propensityExpressions),1);

            for i=1:length(propensityExpressions)
                tmp = propensityExpressions{i};

                if ~isempty(inputs)
                    % Substitute the inputs for their expressions
                    for j = 1:size(inputs,1)
                        nm = inputs{j,1};
                        va = inputs{j,2};
                        tmp = strrep(tmp, nm, ['(' va ')']);
                    end
                end
                
                % Check if propensity contains logical statements.
                logPattern = ["==",">","<"];
                hasLogical=false;
                for j=1:length(logPattern)
                    if contains(tmp,logPattern{j})
                        hasLogical=true;
                    end
                end


                if ~isempty(parDict)
                    parameterCount = length(parDict.keys);
                    parameterNames = parDict.keys;
                    % Substitute the parameters for the values in the propensity
                    % expression
                    if ~hasLogical
                        tmp = str2sym(tmp);
                    end
                    for j = 1:parameterCount
                        nm = parameterNames{j};
                        va = parDict(nm);
                        if ~hasLogical
                            tmp = subs(tmp,nm,num2str(va));
                        else
                            tmp = strrep(tmp, nm, num2str(va));
                        end
                    end
                    tmp = char(tmp);
                end

                if strcmp(str_or_fun,'str')
                    for iSp = 1:length(varNames)
                        tmp = strrep(tmp, varNames{iSp}, ['x(',num2str(iSp),',:)']);
                    end
%                     tmp = strrep(tmp, 'x2', 'x(2,:)');
%                     tmp = strrep(tmp, 'x3', 'x(3,:)');
%                     if ~ismember('x', tmp)
%                         tmp = [tmp '* ones(1, size(x,2))'];
%                     end

                    operators = {'*', '/', '^'};                    
                    for iOp = 1:length(operators)
                        tmp = strrep(tmp, operators{iOp}, ['.' operators{iOp}]);
                    end                      
                    props{i} = tmp;
                elseif strcmp(str_or_fun,'fun')
                    for iSp = 1:length(varNames)
                        tmp = strrep(tmp, varNames{iSp}, ['x(',num2str(iSp),')']);
                    end
%                     tmp = strrep(tmp, 'x1', 'x(1)');
%                     tmp = strrep(tmp, 'x2', 'x(2)');
%                     tmp = strrep(tmp, 'x3', 'x(3)');                    
                    props{i} = str2func(['@(x,t)',tmp]);
                end
            end
        end
    end
end

function dw = findDerivativesSyms(wstr, parNames, varNames)
% Find the derivatives of a propensity function given in wstr with
% respect to the parameters given in parNames. This method first converts
% wstr into a Matlab symbolic expression, then differentiates using Matlab
% Symbolic Toolbox.
% 
% Parameters
% ----------
%
%   wstr: string
%       String expression of the propensity function.
%
%   parNames: cell array of strings
%       Parameter names.
%
% Returns
% -------
%
%   dw: cell array of strings
%     Expressions for the derivatives of the function expressed in `wstr`
%     with respect to the parameters listed in `parNames`.
arguments
    wstr
    parNames
    varNames = {'x1','x2','x3'}
end
num_species = length(varNames);
% Make assumptions on the parameters and variables, such that they are real
num_pars = length(parNames);
for i = 1:num_species
    eval(['syms ',varNames{i}]);
    eval(['assume(',varNames{i},', ''int'')']);
end
for i = 1:length(parNames)
    eval(['syms ',parNames{i}]);
    eval(['assume(',parNames{i},', ''real'')']);
end
syms t;
assume(t, 'real');% assumption for the time variable
% Compute the derivative symbolically
f = str2sym(wstr);
dw = cell(num_pars, 1);
for i = 1:num_pars
    dw{i} = char(diff(f, parNames{i}));
    if string(dw{i}) == "0" || string(dw{i}) == "0.0"
        dw{i} = [];
    end
end
end

function [dfk, ftx, isFactorizable] = separateExpressionParVxt(expr,parameterNames, varNames)
arguments
    expr
    parameterNames
    varNames = {'x1';'x2';'x3'}
end
np = length(parameterNames);
isFactorizable = zeros(1,np);
dfk = cell(1,np);
ftx = cell(1,np);
for ip = 1:np    
    dfk{ip} = [];
    ftx{ip} = [];
    isFactorizable(ip) = false;
    
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
        % Is it a function of only k or of only {x,t}?
        if (~contains(expr, parameterNames{ip}))
            isFactorizable(ip) = true;
            dfk{ip} = [];
            ftx{ip} = 1;
        elseif (~contains(expr, 'x')&&~contains(expr, 't'))
            isFactorizable(ip) = true;
            syms(varNames,'positive')
            syms(parameterNames{ip},'positive')
            eval(['fstr = diff(str2sym(expr),',parameterNames{ip},');']);
            dfk{ip} = (fstr);
            ftx{ip} = 1;
        else
            try
                syms t real
                syms(varNames,'positive')
                syms(parameterNames{ip},'positive')
                eval(['f = str2sym(expr)/',parameterNames{ip},';']);
                fstr = char(f);
                if (~contains(fstr, parameterNames{ip}))
                    isFactorizable(ip) = true;
                    dfk{ip} = 1;
                    ftx{ip} = f;
%                 elseif (~contains(expr, 'x')&&~contains(expr, 't'))
%                     isFactorizable(ip) = true;
%                     fk{ip} = fstr;
%                     fxt{ip} = '1';
                else
                    error('not linear')
                end
            catch
                isFactorizable(ip) = false;
            end
        end
    end
    
%     isFactorizable(ip) = true;
%     f1 = expr(left_bracket_loc(1)+1:right_bracket_loc(1)-1);
%     f2 = expr(left_bracket_loc(2)+1:right_bracket_loc(2)-1);
% 
%     if (contains(f1, 't') || contains(f1, 'I'))
%         ft = f1;
%         fx = f2;
%     else
%         ft = f2;
%         fx = f1;
%     end
end
end
