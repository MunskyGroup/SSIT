%% Spatial Compartement FIM Comparison
clear
close all 
addpath(genpath('../../../src'));

%% Compare Determinate of the Different Models Parameter Uncertanty
% Params
rng('shuffle');
numSamplesToGenerate = 1000;

% Generate Sample for Single Spatial Model
GaussianSampler_4D = @SampleMultiVariableGaussian;
mean = [100, 100 ,100, 100];
var = [50, 50, 50, 50];
GaussianSampler_4D(1, mean, var);
model1Generator = @GenerateSingleSpatialStateModel;
SingleSpatial_Results = GenerateTableFromSampling(model1Generator, GaussianSampler_4D, numSamplesToGenerate, nan);

% Generate Samples for Two Component Spatial Model Pair to the first
GaussianSampler_1D = @SampleMultiVariableGaussian;
mean =  100;
var = 50;
GaussianSampler_1D(size(SingleSpatial_Results,1), mean, var);
model2Generator = @GenerateTwoSpatialStateModel;
TwoSpatial_Results = GenerateTableFromPairSamples(model2Generator, GaussianSampler_1D, [1, 2, 3, 5], SingleSpatial_Results, {'kon', 'koff', 'kr', 'kd'}, 4);

% Generate Samples for Two Component Spatial Model Pair to the first
GaussianSampler_1D = @SampleMultiVariableGaussian;
mean =  100;
var = 50;
GaussianSampler_1D(size(SingleSpatial_Results,1), mean, var);
model3Generator = @GenerateThreeSpatialStateModel;
ThreeSpatial_Results = GenerateTableFromPairSamples(model3Generator, GaussianSampler_1D, [1, 2, 3, 6], TwoSpatial_Results, {'kon', 'koff', 'kr', 'D', 'kd'}, 4);


% Merge Results 
results = join(SingleSpatial_Results, TwoSpatial_Results, Keys={'kon', 'koff', 'kr', 'kd'});
results = join(results, ThreeSpatial_Results, Keys={'kon', 'koff', 'kr', 'kd', 'D'});


%% Parameter Generators
function points = SampleMultiVariableGaussian(varargin)
    persistent numSamples means covMatrix;
    % If this is the first call, initialize the persistent variables
        % If new parameters are provided, update the persistent variables
    if nargin == 3
        numSamples = varargin{1};  % Update number of samples
        means = varargin{2};   % Update mean
        vars = varargin{3};      % Update standard deviation
        covMatrix = diag(vars.^2);  % Diagonal matrix with variances (std^2) on the diagonal
    end

    if ~isempty(numSamples) && ~isempty(covMatrix) && ~isempty(means)
        points = mvnrnd(means, covMatrix, numSamples);
        points(points < 0) = 0.0001; 
    end
end

function points = SampleUniformDistribution(varargin)
    persistent numSamples mins maxs;
    % If this is the first call, initialize the persistent variables
        % If new parameters are provided, update the persistent variables
    if nargin == 3
        numSamples = varargin{1};  % Update number of samples
        mins = varargin{2};
        maxs = varargin{3};
    end

    if ~isempty(numSamples) && ~isempty(mins) && ~isempty(maxs)
        points = rand(numSamples, length(maxs));
        points = mins + (maxs-mins) .* points;
    end
end


%% Model Creation
function model = GenerateSingleSpatialStateModel(varargin)
%{
          kon        kr     kd
gene_off <-> gene_on -> RNA -> 0
          koff
%}
    varargin = varargin{1}; % IDK where the extra cell comes from please explain
    model = SSIT();
    model.parameters = {'kon',varargin(1); 'koff',varargin(2); 'kr',varargin(3); 'kd',varargin(4)};
    model.species = {'gene_on';'gene_off'; 'RNA'};
    model.stoichiometry = [1, -1, 0, 0;...
                        -1, 1, 0, 0;...
                        0, 0, 1, -1];
    model.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'kd*RNA'};
    model.initialCondition = [0;1;0;];
    model.summarizeModel
    model = model.formPropensitiesGeneral('OneComparment_SpatialModel');
end

function model = GenerateTwoSpatialStateModel(varargin)
%{
          kon        kr     D               kd
gene_off <-> gene_on -> RNA -> Distance RNA -> 0
          koff
%}

    model = SSIT();
    varargin = varargin{1}; % IDK where the extra cell comes from please explain
    model.parameters = {'kon',varargin(1);'koff',varargin(2);'kr',varargin(3); 'D', varargin(4); 'kd', varargin(5)};
    model.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'};
    model.stoichiometry = [1, -1, 0, 0, 0;... % species x reaction
                        -1, 1, 0, 0, 0;...
                        0, 0, 1, -1, 0;...
                        0, 0, 0, 1, -1;];
    model.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'D*RNA'; 'kd*D_RNA'};
    model.initialCondition = [0;1;0;0;];
    model.summarizeModel
    model = model.formPropensitiesGeneral('TwoComparment_SpatialModel');
end

function model = GenerateThreeSpatialStateModel(varargin)
%{
          kon        kr              kt         D          kd
gene_off <-> gene_on -> Nascent RNA -> Near RNA -> Far RNA -> 0
          koff
%}

    model = SSIT();
    varargin = varargin{1}; % IDK where the extra cell comes from please explain
    model.parameters = {'kon',varargin(1);'koff',varargin(2);'kr',varargin(3); 'kt',varargin(4); 'D', varargin(5); 'kd', varargin(6)};
    model.species = {'gene_on';'gene_off'; 'nascent_RNA'; 'close_RNA'; 'distant_RNA'};
    model.stoichiometry = [1, -1, 0, 0, 0, 0;...
                           -1, 1, 0, 0, 0, 0;...
                           0, 0, 1, -1, 0, 0;...
                           0, 0, 0, 1, -1, 0;...
                           0, 0, 0, 0, 1, -1;];
    model.propensityFunctions = {'kon*gene_off';'koff*gene_on'; 'kr*gene_on'; 'kt*nascent_RNA'; 'D*close_RNA'; 'kd*distant_RNA'};
    model.initialCondition = [0;1;0;0;0;];
    model.summarizeModel
    model = model.formPropensitiesGeneral('ThreeComparment_SpatialModel');
end


%% FSP Solver, Soles Sensitivity, Computes FIM
function [FIMTotal, FSPsoln, bounds, mleCovEstimate, fimMetrics] = GenerateFSP(model)
    % (2) Solve FSP for model
    model.solutionScheme = 'FSP';    % Set solution scheme to FSP.
    [FSPsoln,model.fspOptions.bounds] = model.solve;  % Solve the FSP analysis

    % (3) Solve FSP Sensitivity
    model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
    [sensSoln,bounds] = model.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem
    
    % (4) Compute FIM using FSP Sensitivity Results
    fimResults = model.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
    cellCounts = 10*ones(size(model.tSpan));  % Number of cells in each experiment.
    [FIMTotal,mleCovEstimate,fimMetrics] = model.evaluateExperiment(fimResults,cellCounts);
end


%% Compute Determinates and making them samples
function determinates = CalcDeterminate(FIMTotal, fimMetrics, J)
    % TODO: Add prior to known parameter
    % TODO: play with fimMetric IDK what these are yet
    % TODO: Return a determinate for each parameter unknown known and
    % combos    
    A = cell2mat(FIMTotal);
    C = inv(A);
     if isnan(J)
        % If J is NaN, use the entire matrix
        B = A;
        cov_unknown = C;
    else
        % Otherwise, use the specified submatrix
        B = A(J, J);
        cov_unknown = C(J, J);
     end
     cov_known = inv(B);

    known_determinate = det(cov_known);
    unknown_determinate = det(cov_unknown);
    determinates = table(known_determinate, unknown_determinate);
end

function sample = GenerateSample(model, J)
    [FIMTotal, FSPsoln, bounds, mleCovEstimate, fimMetrics] = GenerateFSP(model);

    determinates = CalcDeterminate(FIMTotal, fimMetrics, J);

    sample = table(model.parameters{:, 2}, ...
        VariableNames={model.parameters{:, 1}}); % I do not understand matlab indexing {} vs ()
    sample = [sample, determinates];
end


%% Calculate many samples and concat them
function results = GenerateTableFromSampling(modelGenerator, ParameterGenerator, numSampleToGen, J)
    for i = 1:numSampleToGen
        params = ParameterGenerator(); % Generates an array of parameter
        model = modelGenerator(params); % Generates SSIT model
        sample = GenerateSample(model, J);

        if ~exist('results', 'var') % Initial Condition/Preallocate
            width = size(sample, 2);
            results = table(size=[numSampleToGen, width],  VariableNames=sample.Properties.VariableNames, ...
                VariableTypes=table2cell(varfun(@class, sample)));
        end
        results(i,:) = sample;
    end
end

function results = GenerateTableFromPairSamples(modelGenerator, ParameterGenerator, J, t, fixedParams, nonFixedParamsIndexs)
    Params = t(:, fixedParams).Variables;
    newParams = ParameterGenerator();
    counter = 1;
    for p = nonFixedParamsIndexs
        Params = [Params(:, 1:p-1), newParams(:,counter), Params(:, p:end)];
        counter = counter + 1;
    end

    for i = 1:size(Params, 1)
        paramSet = Params(i, :);
        model = modelGenerator(paramSet); % Generates SSIT model
        sample = GenerateSample(model, J);

        if ~exist('results', 'var') % Initial Condition/Preallocate
            width = size(sample, 2);
            results = table(size=[size(fixedParams, 1), width],  VariableNames=sample.Properties.VariableNames, ...
                VariableTypes=table2cell(varfun(@class, sample)));
        end
        results(i,:) = sample;
    end
end
