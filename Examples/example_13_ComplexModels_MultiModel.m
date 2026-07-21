%% SSIT/Examples/example_13_ComplexModels_MultiModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.4.4: Complex models: Multiple Models with shared data/params
%   * Fit multiple models and data sets with shared parameters using
%   SSITMultiModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select scRNA-seq models to include in multimodel:
  ModelNames = {'Model_DUSP1','Model_RUNX1','Model_BIRC3'};

  for i = 1:3
      % Load and prepare codes for each model
        Models{i} = SSIT(['seqModels/',ModelNames{i}]);

      % Set free parameters
        Models{i}.fittingOptions.modelVarsToFit = 1:9;  
        
      % Assign shared (1,2) and unshared (3,..,9) parameters:
        ParInds{i} = [1,2,(3:9)+7*(i-1)];
  end

  Models{1}.formPropensitiesGeneral;

% Set a constraint on model parameters:
  Constraint = @(x) -var(log10([x(3:9);x(7*1+(3:9));x(7*2+(3:9))]));

% Create and initialize multimodel:
  combinedModel = SSITMultiModel(Models, ParInds, Constraint);
  combinedModel = combinedModel.initializeStateSpaces;

% Fit the multimodel:
  fitOptions = optimset('Display','iter','MaxIter',1000);
  combinedModel = combinedModel.maximizeLikelihood(fitOptions=fitOptions);