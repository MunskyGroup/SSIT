function constraintFunctions = readConstraintsForAdaptiveFsp(app,species,Data)
% Updates and creates the constraint boundary functions from the
% constraints table within the FSP GUI
arguments
    app=[];
    species=[];
    Data=[];
end
if isempty(app)
    app.SSITModel.species=species;
    app.FspConstraintTable.Data=Data;
    app.DataLoadingAndFittingTabOutputs.boundIndex=[];
end

nSpecies = length(app.SSITModel.species);
if nSpecies>=1
    spNames = app.SSITModel.species{1};
    for i = 2:nSpecies
        spNames = [spNames,',',app.SSITModel.species{i}];
    end
else
    spNames='';
end

% Updates bounds with the constraints added to the app
new_index = size(app.FspConstraintTable.Data,1);
if isempty(app.DataLoadingAndFittingTabOutputs.boundIndex)
    for i=1:size(app.FspConstraintTable.Data,1)
        Cons = app.FspConstraintTable.Data{i,1};                     % Assigns a variable to represent the value of functions of the constraints function
        app.FspTabOutputs.fConstraints{i} = str2func(['@(',spNames,')',Cons]);
        app.FspTabOutputs.bounds(i) = (app.FspConstraintTable.Data{i,3});
    end
    [~,app.DataLoadingAndFittingTabOutputs.boundIndex] = size(app.FspTabOutputs.bounds);
elseif new_index >= app.DataLoadingAndFittingTabOutputs.boundIndex % these allow if the index changes of the bounds by increasing or decreasing to prevent an error
    for i=1:size(app.FspConstraintTable.Data,1)
        Cons = app.FspConstraintTable.Data{i,1};                     % Assigns a variable to represent the value of functions of the constraints function
        app.FspTabOutputs.fConstraints{i} = str2func(['@(',spNames,')',Cons]);
        app.FspTabOutputs.bounds(i) = (app.FspConstraintTable.Data{i,3});
    end
    app.DataLoadingAndFittingTabOutputs.boundIndex = new_index;
elseif new_index < app.DataLoadingAndFittingTabOutputs.boundIndex
    nConstr = size(app.FspConstraintTable.Data,1);
    B_fun2 = cell(1,nConstr);
    Bounds2 = zeros(1,nConstr);
    for i=1:nConstr
        Cons = app.FspConstraintTable.Data{i,1};                     % Assigns a variable to represent the value of functions of the constraints function
        B_fun2{i} = str2func(['@(',spNames,')',Cons]);
        Bounds2(i) = (app.FspConstraintTable.Data{i,3});
    end
    app.FspTabOutputs.fConstraints = B_fun2;
    app.FspTabOutputs.bounds = Bounds2;
    app.DataLoadingAndFittingTabOutputs.boundIndex = new_index;
end
constraintFunctions = buildTemporaryConstraintFunction(app.FspConstraintTable.Data(:,1), nSpecies);
end

function str2 = convertFormat(str1,numSpecies)
% Convert a string of the form 'x1 * x2' into 'x(1).*x(2)'
str2 = str1;
for i = numSpecies:-1:1
    str2 = strrep(str2, ['x', num2str(i)], ['x(', num2str(i), ',:)']);
end

operators = {'*', '/'};
for i = 1:length(operators)
    str2 = strrep(str2, operators{i}, ['.', operators{i}]);
end

% Cache-friendly rewrite for repeated log terms.
str2 = regexprep(str2, 'log\(\s*x\((\d+)\s*,:\)\s*\+\s*1\s*\)', 'LOGXP1($1,:)');
end

function constraintFunctions = buildTemporaryConstraintFunction(constraintExpressions, numSpecies)
% Create a cached temporary MATLAB function that evaluates all constraints.
cacheKey = strjoin([constraintExpressions; {num2str(numSpecies)}], '||');
persistent functionCache
if isempty(functionCache)
    functionCache = containers.Map('KeyType','char','ValueType','any');
end
if isKey(functionCache, cacheKey)
    constraintFunctions = functionCache(cacheKey);
    return
end

tempName = tempname;
[~, functionName] = fileparts(tempName);
functionName = ['ssit_generated_constraints_', functionName];
functionPath = fullfile(tempdir, [functionName, '.m']);
fileId = fopen(functionPath, 'w');
if fileId < 0
    error('Unable to create temporary constraint function file.')
end

fprintf(fileId, 'function Y = %s(x, b)\n', functionName);
fprintf(fileId, 'LOGXP1 = log(max(0.001,x-3));\n');
% fprintf(fileId, 'LOGXP1 = log(x+1);\n');
if isempty(constraintExpressions)
    fprintf(fileId, 'if nargin < 2\n');
    fprintf(fileId, 'Y = zeros(0, size(x,2));\n');
    fprintf(fileId, 'else\n');
    fprintf(fileId, 'Y = false(0, size(x,2));\n');
    fprintf(fileId, 'end\n');
else
    formattedExpressions = cell(size(constraintExpressions));
    formattedMaskExpressions = cell(size(constraintExpressions));
    for i = 1:numel(constraintExpressions)
        formattedExpressions{i} = convertFormat(constraintExpressions{i}, numSpecies);
        formattedMaskExpressions{i} = ['(', formattedExpressions{i}, ') > b(', num2str(i), ')'];
    end
    fprintf(fileId, 'if nargin < 2\n');
    fprintf(fileId, 'Y = [%s];\n', strjoin(formattedExpressions, '; '));
    fprintf(fileId, 'else\n');
    fprintf(fileId, 'b = b(:);\n');
    fprintf(fileId, 'Y = [%s];\n', strjoin(formattedMaskExpressions, '; '));
    fprintf(fileId, 'end\n');
end
fprintf(fileId, 'end\n');
fclose(fileId);

if ~contains([path, pathsep], [tempdir, pathsep])
    addpath(tempdir);
end
constraintFunctions = str2func(functionName);
functionCache(cacheKey) = constraintFunctions;
end