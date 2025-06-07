function configurations = multiplyConfigurations(configurables)
    arguments
        configurables (1, :) AbstractExperimentConfigurable
    end
    
    if isempty(configurables)
        configurations = [];
    else
        configurations = ExperimentConfiguration();
        if isscalar(configurables)
            multipliedConfigurables = configurables.multiply();          
            configurations.Configurables = multipliedConfigurables;
        else
            % We operate recursively by first multiplying all configurables
            % except the last one; suppose that this results in M
            % configurations, whereas the last configurable has N possible
            % values. Then we will have M*N configurations after
            % multiplying by the last configurable.

            % In the recursive step, we need to extract the Configurables
            % property (a list of AbstractExperimentConfigurable) because
            % this method always returns an ExperimentConfiguration.

            preMultipliedConfigurables = ...
                multiplyConfigurations(configurables(1:end-1));          
            preMultipliedConfigurables = ...
                preMultipliedConfigurables.Configurables;
            M = length(preMultipliedConfigurables);

            lastConfigurables = configurables(end).multiply();
            N = length(lastConfigurables);           

            configurations = ...
                createArray(1, M * N, "ExperimentConfiguration");
            for configurableIdx = 0:(M * N - 1)
                % We will create a "grid" in which the row corresponds to
                % the list of configurables from the "pre-multiplied" ones
                % and the column to the "last" configurable. The math for
                % doing this is straightforward when using zero-based
                % indexing, so we will use that and simply increment when 
                % indexing into the configurables themselves.

                rowIdx = idivide(configurableIdx, int32(M));
                colIdx = mod(configurableIdx, N);
                
                configurables = ...
                    [preMultipliedConfigurables(rowIdx + 1) ...
                    lastConfigurables(colIdx + 1)];
                
                configurations(configurableIdx + 1).Configurables = ...
                    configurables;
            end            
        end % [configurables contains 2+]
    end % [configurables is non-empty]
end