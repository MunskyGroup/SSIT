function constraintFunctions = readConstraintsForAdaptiveFsp(app,species,Data)
% Updates and creates the constraint boundary functions from the
% constraints table within the FSP GUI
arguments
    app=[];
    species=[];
    Data=[];
end
if isempty(app)
    app.ReactionsTabOutputs.varNames=species;
    app.FspConstraintTable.Data=Data;
    app.DataLoadingAndFittingTabOutputs.boundIndex=[];
end   


nSpecies = length(app.ReactionsTabOutputs.varNames);
spNames = app.ReactionsTabOutputs.varNames{1};
for i = 2:nSpecies
    spNames = [spNames,',',app.ReactionsTabOutputs.varNames{i}];
end

% Updates bounds with the constraints added to the app
new_index = size(app.FspConstraintTable.Data,1);
fstr = '[';
if isempty(app.DataLoadingAndFittingTabOutputs.boundIndex)
    for i=1:size(app.FspConstraintTable.Data,1)
        Cons = app.FspConstraintTable.Data{i,1};                     % Assigns a variable to represent the value of functions of the constraints function        
        app.FspTabOutputs.fConstraints{i} = str2func(['@(',spNames,')',Cons]);
        app.FspTabOutputs.bounds(i) = (app.FspConstraintTable.Data{i,3});
        fstr = [ fstr ' ' Cons ';'];
    end
    [~,app.DataLoadingAndFittingTabOutputs.boundIndex] = size(app.FspTabOutputs.bounds);
    fstr = [fstr ' ]'];
elseif new_index >= app.DataLoadingAndFittingTabOutputs.boundIndex % these allow if the index changes of the bounds by increasing or decreasing to prevent an error
    for i=1:size(app.FspConstraintTable.Data,1)
        Cons = app.FspConstraintTable.Data{i,1};                     % Assigns a variable to represent the value of functions of the constraints function
        app.FspTabOutputs.fConstraints{i} = str2func(['@(',spNames,')',Cons]);
        app.FspTabOutputs.bounds(i) = (app.FspConstraintTable.Data{i,3});
        fstr = [ fstr ' ' Cons ';'];
    end
    app.DataLoadingAndFittingTabOutputs.boundIndex = new_index;
    fstr = [fstr ' ]'];
elseif new_index < app.DataLoadingAndFittingTabOutputs.boundIndex
    for i=1:size(app.FspConstraintTable.Data,1)
        Cons = app.FspConstraintTable.Data{i,1};                     % Assigns a variable to represent the value of functions of the constraints function
        B_fun2{i} = str2func(['@(',spNames,')',Cons]);
        Bounds2(i) = (app.FspConstraintTable.Data{i,3});
        fstr = [ fstr ' ' Cons ';'];
    end
    app.FspTabOutputs.fConstraints = B_fun2;
    app.FspTabOutputs.bounds = Bounds2;
    app.DataLoadingAndFittingTabOutputs.boundIndex = new_index;
    fstr = [fstr ' ]'];
end
fstr = convertFormat(fstr,nSpecies);
constraintFunctions = str2func(['@(x) ' fstr]);
end

function str2 = convertFormat(str1,numSpecies)
% Convert a string of the form 'x1 * x2' into 'x(1).*x(2)'
str2 = str1;
for i = 1:numSpecies
    str2 = strrep(str2, ['x', num2str(i)], ['x(', num2str(i), ',:)']);
end

operators = {'*', '/'};
for i = 1:length(operators)
    str2 = strrep(str2, operators{i}, ['.', operators{i}]);
end
end