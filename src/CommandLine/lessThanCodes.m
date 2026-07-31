function result = lt(obj1, obj2)
arguments
    obj1 (1, 1) MBExperimentInputConfigurable
    obj2 (1, 1) MBExperimentInputConfigurable
end

if obj1.InputName < obj2.InputName
    result = true;
elseif obj1.InputName > obj2.InputName
    result = false;
else
    % The input names are the same. The smaller one will be the
    % one that has fewer values. If each object has the same
    % number of values, iterate over each pair of elements to
    % determine the smaller.
    if obj1.numberOfValues() < obj2.numberOfValues()
        result = true;
    elseif obj1.numberOfValues() > obj2.numberOfValues()
        result = false;
    else
        result = false;
        for valIdx = 1:obj1.numberOfValues()
            if obj1.Values(valIdx) < obj2.Values(valIdx)
                result = true;
                break;
            elseif obj1.Values(valIdx) > obj2.Values(valIdx)
                result = false;
                break;
            end
        end
    end
end
end % lt MBExperimentInputConfigurable

function result = lt(obj1, obj2)
arguments
    obj1 (1, 1) MBExperimentParameterConfigurable
    obj2 (1, 1) MBExperimentParameterConfigurable
end

if obj1.ParameterName < obj2.InputName
    result = true;
elseif obj1.ParameterName > obj2.InputName
    result = false;
else
    % The input names are the same. The smaller one will be the
    % one that has fewer values. If each object has the same
    % number of values, iterate over each pair of elements to
    % determine the smaller.
    if obj1.numberOfValues() < obj2.numberOfValues()
        result = true;
    elseif obj1.numberOfValues() > obj2.numberOfValues()
        result = false;
    else
        result = false;
        for valIdx = 1:obj1.numberOfValues()
            if obj1.Values(valIdx) < obj2.Values(valIdx)
                result = true;
                break;
            elseif obj1.Values(valIdx) > obj2.Values(valIdx)
                result = false;
                break;
            end
        end
    end
end
end % lt MBExperimentParameterConfigurable