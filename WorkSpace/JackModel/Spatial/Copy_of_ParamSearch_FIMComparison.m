%% Spatial Compartement FIM Comparison
clear
close all 
addpath(genpath('../../../src'));

%% Compare Determinate of the Different Models Parameter Uncertanty
% Params
rng('shuffle');
numSamplesToGenerate = 100;

% Generate Sample for Single Spatial Model
GaussianSampler_4D = @SampleMultiVariableGaussian;
mean = [10, 1 ,1, 0.1];
var = [3, 0.25, 0.25, 0.01];
GaussianSampler_4D(1, mean, var);
model1Generator = @GenerateSingleSpatialStateModel;
SingleSpatial_Results = GenerateTableFromSampling(model1Generator, GaussianSampler_4D, numSamplesToGenerate, nan);

% Generate Samples for Two Component Spatial Model Pair to the first
GaussianSampler_1D = @SampleMultiVariableGaussian;
mean =  1;
var = 0.25;
GaussianSampler_1D(size(SingleSpatial_Results,1), mean, var);
model2Generator = @GenerateTwoSpatialStateModel;
TwoSpatial_Results = GenerateTableFromPairSamples(model2Generator, GaussianSampler_1D, 4, SingleSpatial_Results, {'kon', 'koff', 'kr', 'kd'});

% Generate Samples for Two Component Spatial Model Pair to the first
GaussianSampler_1D = @SampleMultiVariableGaussian;
mean =  1;
var = 0.25;
GaussianSampler_1D(size(SingleSpatial_Results,1), mean, var);
model3Generator = @GenerateThreeSpatialStateModel;
ThreeSpatial_Results = GenerateTableFromPairSamples(model3Generator, GaussianSampler_1D, 4, SingleSpatial_Results, {'kon', 'koff', 'kr', 'kd'});


% Merge Results 
results = join(SingleSpatial_Results, TwoSpatial_Results, Keys={'kon', 'koff', 'kr', 'kd'});


%% Display Results
diffusion2decayRatio = results.D ./ results.kd;

fig = figure();
scatter(results.known_param_determinate_SingleSpatial_Results, results.known_param_determinate_TwoSpatial_Results, ...
        100, diffusion2decayRatio, 'filled');
colorbar;
xlabel('|E| w/o spatial');
ylabel('|E| w/ spatial');
title('Determinate of Parameter Confidence with Known Diffusion Coeffient');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
xRange = linspace(min(results.known_param_determinate_SingleSpatial_Results), max(results.known_param_determinate_SingleSpatial_Results), 100);
plot(xRange, xRange, 'r-', 'LineWidth', 2);
hold off;
saveas(fig, 'KnownDiffusion_Determinates.png');


fig = figure();
scatter(results.unknown_param_determinate_SingleSpatial_Results, results.unknown_param_determinate_TwoSpatial_Results, ...
        100, diffusion2decayRatio, 'filled');
colorbar;
xlabel('|E| w/o spatial');
ylabel('|E| w/ spatial');
title('Determinate of Parameter Confidence with Unknown Diffusion Coeffient');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold on;
xRange = linspace(min(results.unknown_param_determinate_SingleSpatial_Results), max(results.unknown_param_determinate_SingleSpatial_Results), 100);
plot(xRange, xRange, 'r-', 'LineWidth', 2);
hold off;
saveas(fig, 'UnknownDiffusion_Determinates.png');




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
    model = model.formPropensitiesGeneral('1Comparment_SpatialModel');
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
    model = model.formPropensitiesGeneral('2Comparment_SpatialModel');
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
    model = model.formPropensitiesGeneral('3Comparment_SpatialModel');
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
function [known_param_determinate, unknown_param_determinate] = CalcDeterminate(FIMTotal, fimMetrics, paramOfInterestIndex)
    % TODO: Add prior to known parameter
    % TODO: play with fimMetric IDK what these are yet
    % TODO: Return a determinate for each parameter unknown known and
    % combos
    A = cell2mat(FIMTotal);
    j = paramOfInterestIndex;
    
    % Known and Unknown Parameter determinate
    covar_A = inv(A);
    if isnan(paramOfInterestIndex)
        B = A;
        C = covar_A;
    else
        B = A([1:j-1, j+1:end], [1:j-1, j+1:end]);
        C = covar_A([1:j-1, j+1:end], [1:j-1, j+1:end]);
    end
    covar_B = inv(B);

    known_param_determinate = det(covar_B);
    unknown_param_determinate = det(C);
end

function sample = GenerateSample(model, paramOfInterestIndex)
    [FIMTotal, FSPsoln, bounds, mleCovEstimate, fimMetrics] = GenerateFSP(model);

    [known_param_determinate, unknown_param_determinate] = CalcDeterminate(FIMTotal, fimMetrics, paramOfInterestIndex);

    sample = table(model.parameters{:, 2}, known_param_determinate, unknown_param_determinate, ...
        VariableNames={model.parameters{:, 1}, 'known_param_determinate', 'unknown_param_determinate'}); % I do not understand matlab indexing {} vs ()
end


%% Calculate many samples and concat them
function results = GenerateTableFromSampling(modelGenerator, ParameterGenerator, numSampleToGen, paramOfInterestIndex)
    for i = 1:numSampleToGen
        params = ParameterGenerator(); % Generates an array of parameter
        model = modelGenerator(params); % Generates SSIT model
        sample = GenerateSample(model, paramOfInterestIndex);

        if ~exist('results', 'var') % Initial Condition/Preallocate
            width = size(sample, 2);
            results = table(size=[numSampleToGen, width],  VariableNames=sample.Properties.VariableNames, ...
                VariableTypes=table2cell(varfun(@class, sample)));
        end

        results(i,:) = sample;
    end
end

function results = GenerateTableFromPairSamples(modelGenerator, ParameterGenerator, paramOfInterestIndex, t, params)
    params = t(:, params).Variables;
    paramOfInterest = ParameterGenerator();
    % TODO: Make it work for mutiple paramOfInterestIndex
    params = [params(:, 1:paramOfInterestIndex-1), paramOfInterest, params(:, paramOfInterestIndex:end)];

    for i = 1:size(params, 1)
        paramSet = params(i, :);
        model = modelGenerator(paramSet); % Generates SSIT model
        sample = GenerateSample(model, paramOfInterestIndex);

        if ~exist('results', 'var') % Initial Condition/Preallocate
            width = size(sample, 2);
            results = table(size=[size(params, 1), width],  VariableNames=sample.Properties.VariableNames, ...
                VariableTypes=table2cell(varfun(@class, sample)));
        end
        results(i,:) = sample;
    end
end




















