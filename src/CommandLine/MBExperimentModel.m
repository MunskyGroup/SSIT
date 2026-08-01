classdef MBExperimentModel < SSIT
    %MBExperimentModel encapsulates a model on which sequential experiment
    %design can be based. It extends the SSIT class with certain properties
    %needed for experiment design.

    properties
        Property1
    end

    methods
        function isEqual = eq(obj1, obj2)
            %This method implements an equality comparator for the
            %MBExperimentModel class. It validates that certain properties
            %match between the two models. It does NOT require exact
            %identity between the models, because this descendant class is
            %only intended for experiment design, and certain properties
            %are either overwritten in the design process (e.g., fitting
            %options) or are unsupported.

            arguments
                obj1 (1, 1) MBExperimentModel
                obj2 (1, 1) MBExperimentModel                
            end

            % We assume equality until proven otherwise and use
            % short-circuiting to rapidly exit if and when inequality is
            % shown.

            isEqual = true; 

            % Parameters must match in name and value:

            isEqual = isEqual && ...
                size(obj1.parameters) == size(obj2.parameters);
            isEqual = isEqual && all(...
                strcmp(obj1.parameters{:, 1}, obj2.parameters{:, 1}));
            isEqual = isEqual && all(...
                obj1.parameters{:, 2} == obj2.parameters{:, 2});

            % Species must match in name:

            isEqual = isEqual && ...
                (size(obj1.species) == size(obj2.species));
            isEqual = isEqual && all(strcmp(obj1.species, obj2.species));

            % Stoichiometries must be identical element-wise:

            isEqual = isEqual && ...
                size(obj1.stoichiometry) == size(obj2.stoichiometry);
            isEqual = isEqual && all(...
                obj1.stoichiometry == obj2.stoichiometry);

            % Propensity functions must match in value.
            % The corresponding executable functions (cell arrays of 
            % SSIT.Propensity) are presumed to match as long as the
            % propensity functions match, so the following properties are
            % ignored:
            %
            % propensitiesGeneral
            % propensitiesGeneralODE
            % propensitiesGeneralMean
            % propensitiesGeneralMeanJac
            % propensitiesGeneralMoments
            % propensitiesGeneralMomentsJac
            % propensityFilePrefix

            isEqual = isEqual && ...
                (size(obj1.propensityFunctions) == ...
                size(obj2.propensityFunctions));
            isEqual = isEqual && all(strcmp(...
                obj1.propensityFunctions, obj2.propensityFunctions));

            % Input expressions must match in name and value:

            isEqual = isEqual && ...
                size(obj1.inputExpressions) == size(obj2.inputExpressions);
            isEqual = isEqual && all(strcmp(...
                obj1.inputExpressions{:, 1}, ...
                obj2.inputExpressions{:, 1}));
            isEqual = isEqual && all(strcmp(...
                obj1.inputExpressions{:, 2}, ...
                obj2.inputExpressions{:, 2}));

            % Constraint functions for FSP must match in value:

            isEqual = isEqual && ...
                (size(obj1.customConstraintFuns) == ...
                size(obj2.customConstraintFuns));
            isEqual = isEqual && all(strcmp(...
                obj1.customConstraintFuns, obj2.customConstraintFuns));

            % Initial population sizes must be identical element-wise:

            isEqual = isEqual && ...
                size(obj1.initialCondition) == size(obj2.initialCondition);
            isEqual = isEqual && all(...
                obj1.initialCondition == obj2.initialCondition);

            % Initial state probability masses must be identical 
            % element-wise:

            isEqual = isEqual && ...
                size(obj1.initialProbs) == size(obj2.initialProbs);
            isEqual = isEqual && all(...
                obj1.initialProbs == obj2.initialProbs);

            % Initial times must be identical:

            isEqual = isEqual && (obj1.initialTime == obj2.initialTime);

            % The following are currently UNSUPPORTED in the design process
            % and therefore must be empty/false:
            %
            % specialEvents
            % PDO
            % Model reduction

            isEqual = isEqual && ...
                isempty(obj1.specialEvents) && isempty(obj2.specialEvents);

            isEqual = isEqual && ...
                isempty(obj1.pdoOptions.PDO) && ...
                isempty(obj2.pdoOptions.PDO) && ...
                isempty(obj1.pdoOptions.unobservedSpecies) && ...
                isempty(obj2.pdoOptions.unobservedSpecies);

            isEqual = isEqual && ...
                ~obj1.modelReductionOptions.useModReduction && ...
                ~obj2.modelReductionOptions.useModReduction;

            % The following are OVERWRITTEN in the design process as
            % necessary and therefore ignored in equality comparison:
            %
            % fspOptions
            % sensOptions
            % ssaOptions
            % fittingOptions
            % tSpan
            % solutionScheme
            % odeIntegrator
            % solutionSchemes
            % dataSet 
            % Solutions

            % The following are IRRELEVANT to the design process and
            % therefore ignored in equality comparison:
            %
            % description
            % GUIProps
               
        % Option to use hybrid model (deterministic + stochastic species)
        %   default: false
        useHybrid = false
        % Struct to define which species of the hybrid model will be
        % modelled using ODEs
        %   default: struct('upstreamODEs',[]);
        %   example: Model.hybridOptions.upstreamODEs = {'offGene','onGene'};
        hybridOptions = struct('upstreamODEs',[]);
        
        end % eq

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end

        function obj = set.Property1(obj, value)
            obj.fittingOptions.logPriorCovariance = value;
        end
    end
end