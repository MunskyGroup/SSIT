function d = createLinkedSpeciesDictionary(...
    species, model)
%createLinkedSpeciesDictionary converts an n-by-2 or n-by-3 cell array
%representing linkages between model and data species into an equivalent
%dictionary mapping model species (keys) to data species (values). If a 
%model is provided, an error will be raised if any model species in the 
%input array is not a model species of the input model.
arguments (Input)
    species (:, :) cell {mustBeValidLinkedSpeciesArray}
    model SSIT {mustBeScalarOrEmpty} = createArray(0, "SSIT")
end

d = configureDictionary("string", "string");

[numberOfSpecies, ~] = size(species);
for speciesIdx = 1:numberOfSpecies
    curModelSpecies = string(species{speciesIdx, 1});
    if ~isempty(model) & ~any(contains(model.species, curModelSpecies))
        error("Could not find linked species " + curModelSpecies + ...
            " in the existing model species: " + ...
            join(string(model.species), ", "))
    end
    curMeasuredDataSpecies = string(species{speciesIdx, 2});
    curComputedDataSpecies = string(species{speciesIdx, 3});

    % The mustBeValidLinkedSpeciesArray method validates that for each
    % model species, exclusively either a measured data species or a
    % computed data species will be supplied.

    if strlength(curMeasuredDataSpecies) > 0
        d(curModelSpecies) = curMeasuredDataSpecies;
    else
        d(curModelSpecies) = curComputedDataSpecies;
    end
end

end

function mustBeValidLinkedSpeciesArray(array)
    mustBeMatrix(array); % Ensure that the array is two-dimensional
    [numberOfRows, numberOfColumns] = size(array);

    if numberOfColumns < 2
        error("A linked species array contains too few columns: " + ...
            numberOfColumns)
    elseif numberOfColumns > 3
        error("A linked species array contains too many columns: " + ...
            numberOfColumns)
    end

    % Validate that every model species name (column 1) is unique:
    if length(unique(array{:, 1})) ~= numberOfRows
        error("A linked species array contains duplicate model species")
    end

    % Validate that every model species has a corresponding entry in EITHER
    % the measured column (column 2) OR the computed column (column 3):

    for speciesIdx = 1:numberOfRows
        measuredVariable = string(array{speciesIdx, 2});
        computedVariable = string(array{speciesIdx, 2});
        hasMeasuredVariable = strlength(measuredVariable) > 0;
        hasComputedVariable = strlength(computedVariable) > 0;
        if hasMeasuredVariable && hasComputedVariable
            error("A linked species array has multiple mappings for " + ...
                "model species " + string(array{speciesIdx, 1}) + ...
                ": " + measuredVariable + " and " + computedVariable)
        elseif ~hasMeasuredVariable && ~hasComputedVariable
            error("A linked species array has no mappings for " + ...
                "model species " + string(array{speciesIdx, 1}))
        end
    end    
end