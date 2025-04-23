%% Generates D vs |FIM| plot
clear
close all 
addpath(genpath('../../src'));
%% Pretext 
% Iyer-Biswas, Srividya, et al. “Stochasticity of Gene Products from Transcriptional Pulsing.”
% This paper explored the parameter space of the one compartment model,
% some key finding was that youu can fix everything relative to one
% constant (kd) and explore the parameters, kr/kd conrols the extent of the
% distribution and not the shape, kon and koff primarly control the shape
% we will explore the five regimes:
% regime 1 - kd = kon = koff
% regime 2 - kon, koff > kd
% regime 3 - kon, koff < kd
% regime 4 - kon < kd < koff
% regime 5 - koff < kd < kon


%% Build Model
% One Compartment Model
OneCompartmentModel = SSIT();
OneCompartmentModel.parameters = {'kon',1; 'koff',1; 'kr',100; 'kd',1};
OneCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'};
OneCompartmentModel.stoichiometry = [1, -1, 0, 0;...
                                    -1, 1, 0, 0;...
                                    0, 0, 1, -1];
OneCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'kd*RNA'};
OneCompartmentModel.initialCondition = [0;1;0;];
OneCompartmentModel.summarizeModel
OneCompartmentModel = OneCompartmentModel.formPropensitiesGeneral('OneComparment_SpatialModel');
OneCompartmentModel.tSpan = linspace(0,10,21);
% OneCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on';'gene_off';};
save('OneComparment_SpatialModel','OneCompartmentModel')

% Two Compartment Model
TwoCompartmentModel = SSIT();
TwoCompartmentModel.parameters = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
TwoCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'};
TwoCompartmentModel.stoichiometry = [1, -1, 0, 0, 0, 0, 0;... % species x reaction
                                    -1, 1, 0, 0, 0, 0, 0;...
                                    0, 0, 1, -1, 0, -1, 1;...
                                    0, 0, 0, 1, -1, 0, -1;];
TwoCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'D*RNA'; 'kd*D_RNA'; 'kd*RNA'; 'D*D_RNA';};
TwoCompartmentModel.initialCondition = [0;1;0;0;];
TwoCompartmentModel.summarizeModel
TwoCompartmentModel = TwoCompartmentModel.formPropensitiesGeneral('TwoComparment_SpatialModel');
TwoCompartmentModel.tSpan = linspace(0,10,21);
% TwoCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on';'gene_off';};
save('TwoComparment_SpatialModel','TwoCompartmentModel')


%% Pipeline
run_all = false;

x = zeros([1, 5]);
y1 = zeros([1, 5]);
y2 = zeros([1, 5]);

% specify pipeline parameters
pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = OneCompartmentModel.parameters;
pipelineArgs.makePlots = false;

% run pipeline on model
saveName = 'OneCompartmentModel_regime1';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
x(1) = outputs.determinates{1, "known_determinate"};


% specify pipeline parameters
pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = 4;
pipelineArgs.pars = TwoCompartmentModel.parameters;
pipelineArgs.makePlots = true;

% run pipeline on model
saveName = 'TwoCompartmentModel_regime1';
if run_all
SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
y1(1) = outputs.determinates{1, "known_determinate"};
y2(1) = outputs.determinates{1, "unknown_determinate"};


%% Run Multiple Regimes 
% regime 2 - kon, koff > kd
pipelineArgs.makePlots = false;
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',1.2; 'koff',1.5; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime2';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
x(2) = outputs.determinates{1, "known_determinate"};

pipelineArgs.param_of_interest_index = 4;
pipelineArgs.pars = {'kon',1.2;'koff', 1.5;'kr',100; 'D', 1; 'kd', 1};
saveName = 'TwoCompartmentModel_regime2';
if run_all
SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
y1(2) = outputs.determinates{1, "known_determinate"};
y2(2) = outputs.determinates{1, "unknown_determinate"};

% regime 3 - kon, koff < kd
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',0.5; 'koff',0.5; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime3';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
x(3) = outputs.determinates{1, "known_determinate"};

pipelineArgs.param_of_interest_index = 4;
pipelineArgs.pars = {'kon',0.5;'koff', 0.5;'kr',100; 'D', 1; 'kd', 1};
saveName = 'TwoCompartmentModel_regime3';
if run_all
SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
y1(3) = outputs.determinates{1, "known_determinate"};
y2(3) = outputs.determinates{1, "unknown_determinate"};

% regime 4 - kon < kd < koff
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',0.5; 'koff',2; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime4';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
x(4) = outputs.determinates{1, "known_determinate"};

pipelineArgs.param_of_interest_index = 4;
pipelineArgs.pars = {'kon',0.5;'koff', 2;'kr',100; 'D', 1; 'kd', 1};
saveName = 'TwoCompartmentModel_regime4';
if run_all
SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
y1(4) = outputs.determinates{1, "known_determinate"};
y2(4) = outputs.determinates{1, "unknown_determinate"};

% regime 5 - koff < kd < kon
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',1.5; 'koff',0.5; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime5';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
x(5) = outputs.determinates{1, "known_determinate"};

pipelineArgs.param_of_interest_index = 4;
pipelineArgs.pars = {'kon',1.5;'koff', 0.5;'kr',100; 'D', 1; 'kd', 1};
saveName = 'TwoCompartmentModel_regime5';
if run_all
SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
y1(5) = outputs.determinates{1, "known_determinate"};
y2(5) = outputs.determinates{1, "unknown_determinate"};

% Plot
figure()
scatter(x, y1, 60, 'r', 'filled');
hold on;
scatter(x, y2, 60, 'b', 'filled');

x_ref = [min(x), max(x)];
y_ref = x_ref;
plot(x_ref, y_ref, 'k--', 'LineWidth', 1.5);

set(gca, 'XScale', 'log', 'YScale', 'log');
legend('Diffusion Coeffient is known', 'Diffusion Coeffient is unknown');
xlabel('|CoV| for one compartment model');
ylabel('|CoV| for two compartment model');
title('|CoV| Comparison in Standard Regimes');


%% Run Multiple Diffusion Coeffients
run_all=false;
D = [0, 0.01, 0.1, 1, 10];
f = zeros([2, 5, length(D)]);
g = zeros([5, length(D)]);

% regime 1 - kd = kon = koff
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',1; 'koff',1; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime1';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
g(1,:) = outputs.determinates{1, "known_determinate"};

for d=1:length(D)
    pipelineArgs.param_of_interest_index = 4;
    pipelineArgs.pars = {'kon',1;'koff',1;'kr',100; 'D', D(d); 'kd', 1};
    saveName = ['TwoCompartmentModel_regime1', '_', num2str(fix(D(d)*100))];
    if run_all
    SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
    end
    load([saveName, '.mat'])
    f(1,1,d) = outputs.determinates{1, "known_determinate"};
    f(2,1,d) = outputs.determinates{1, "unknown_determinate"};
end

% regime 2 - kon, koff > kd
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',1.2; 'koff',1.5; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime2';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
g(2,:) = outputs.determinates{1, "known_determinate"};

for d=1:length(D) 
    pipelineArgs.param_of_interest_index = 4;
    pipelineArgs.pars = {'kon',1.2;'koff', 1.5;'kr',100; 'D', D(d); 'kd', 1};
    saveName = ['TwoCompartmentModel_regime2', '_', num2str(fix(D(d)*100))];
    if run_all
    SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
    end
    load([saveName, '.mat'])
    f(1,2,d) = outputs.determinates{1, "known_determinate"};
    f(2,2,d) = outputs.determinates{1, "unknown_determinate"};
end

% regime 3 - kon, koff < kd
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',0.5; 'koff',0.5; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime3';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
g(3,:) = outputs.determinates{1, "known_determinate"};

for d=1:length(D)

    pipelineArgs.param_of_interest_index = 4;
    pipelineArgs.pars = {'kon',0.5;'koff', 0.5;'kr',100; 'D', D(d); 'kd', 1};
    saveName = ['TwoCompartmentModel_regime3', '_', num2str(fix(D(d)*100))];
    if run_all
    SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
    end
    load([saveName, '.mat'])
    f(1,3,d) = outputs.determinates{1, "known_determinate"};
    f(2,3,d) = outputs.determinates{1, "unknown_determinate"};
end

% regime 4 - kon < kd < koff
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',0.5; 'koff',2; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime4';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
g(4,:) = outputs.determinates{1, "known_determinate"};
    
for d=1:length(D)
    pipelineArgs.param_of_interest_index = 4;
    pipelineArgs.pars = {'kon',0.5;'koff', 2;'kr',100; 'D', D(d); 'kd', 1};
    saveName = ['TwoCompartmentModel_regime4', '_', num2str(fix(D(d)*100))];
    if run_all
    SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
    end
    load([saveName, '.mat'])
    f(1,4,d) = outputs.determinates{1, "known_determinate"};
    f(2,4,d) = outputs.determinates{1, "unknown_determinate"};
end

% regime 5 - koff < kd < kon
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = {'kon',1.5; 'koff',0.5; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime5';
if run_all
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
end
load([saveName, '.mat'])
g(5,:) = outputs.determinates{1, "known_determinate"};

for d=1:length(D)
    pipelineArgs.param_of_interest_index = 4;
    pipelineArgs.pars = {'kon',1.5;'koff', 0.5;'kr',100; 'D', D(d); 'kd', 1};
    saveName = ['TwoCompartmentModel_regime5', '_', num2str(fix(D(d)*100))];
    if run_all
    SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
    end
    load([saveName, '.mat'])
    f(1,5,d) = outputs.determinates{1, "known_determinate"};
    f(2,5,d) = outputs.determinates{1, "unknown_determinate"};
end

%% Plot
figure()
scatter(D, squeeze(f(1,1,:)), 60, 'r', 'filled'); % regime 1 known
hold on;
scatter(D, squeeze(f(2,1,:)), 60, 'r'); % regime 1 unknown
plot(D, squeeze(g(1,:)), 'r--', 'LineWidth', 1.5);

scatter(D, squeeze(f(1,2,:)), 60, 'g', 'filled'); % regime 2 known
scatter(D, squeeze(f(2,2,:)), 60, 'g'); % regime 2 unknown
plot(D, squeeze(g(2,:)), 'g--', 'LineWidth', 1.5);

scatter(D, squeeze(f(1,3,:)), 60, 'b', 'filled'); % regime 3 known
scatter(D, squeeze(f(2,3,:)), 60, 'b'); % regime 3 unknown
plot(D, squeeze(g(3,:)), 'b--', 'LineWidth', 1.5);

scatter(D, squeeze(f(1,4,:)), 60, 'c', 'filled'); % regime 4 known
scatter(D, squeeze(f(2,4,:)), 60, 'c'); % regime 4 unknown
plot(D, squeeze(g(4,:)), 'c--', 'LineWidth', 1.5);

scatter(D, squeeze(f(1,5,:)), 60, 'm', 'filled'); % regime 5 known
scatter(D, squeeze(f(2,5,:)), 60, 'm'); % regime 5 unknown
plot(D, squeeze(g(5,:)), 'm--', 'LineWidth', 1.5);

set(gca, 'XScale', 'log', 'YScale', 'log');
legend('Regime 1 - known', 'regime 1 - unknown', 'regime 1 - one compartment', ...
    'Regime 2 - known', 'regime 2 - unknown', 'regime 2 - one compartment', ...
    'Regime 3 - known', 'regime 3 - unknown', 'regime 3 - one compartment', ...
    'Regime 4 - known', 'regime 4 - unknown', 'regime 4 - one compartment', ...
    'Regime 5 - known', 'regime 5 - unknown', 'regime 5 - one compartment' ...
    );
xlabel('D');
ylabel('|CoV|');
title('D vs |CoV|');








