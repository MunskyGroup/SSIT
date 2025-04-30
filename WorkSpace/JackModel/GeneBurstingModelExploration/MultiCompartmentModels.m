%% Generates D vs |FIM| plot
clear
close all 
addpath(genpath('../../../src'));
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

% We will only look at regime 1 for this case and test multiple compartment
% models, with different transition rates


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
OneCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on','gene_off'};
save('OneComparment_SpatialModel','OneCompartmentModel')

% Two Compartment Model
TwoCompartmentModel = SSIT();
TwoCompartmentModel.parameters = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
TwoCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'};
TwoCompartmentModel.stoichiometry = [1, -1, 0, 0, 0, 0, 0;... % species x reaction
                                    -1, 1, 0, 0, 0, 0, 0;...
                                    0, 0, 1, -1, 0, -1, 1;...
                                    0, 0, 0, 1, -1, 0, -1;];
TwoCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'D*RNA'; 'kd*D_RNA'; 'kd*RNA'; 'D*D_RNA'};
TwoCompartmentModel.initialCondition = [0;1;0;0;];
TwoCompartmentModel.summarizeModel
TwoCompartmentModel = TwoCompartmentModel.formPropensitiesGeneral('TwoComparment_SpatialModel');
TwoCompartmentModel.tSpan = linspace(0,10,21);
TwoCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on','gene_off'};
save('TwoComparment_SpatialModel','TwoCompartmentModel')

% Three Compartment Model
ThreeCompartmentModel = SSIT();
ThreeCompartmentModel.parameters = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
ThreeCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'; 'D1_RNA'};
ThreeCompartmentModel.stoichiometry = [1,-1, 0, 0, 0, 0, 0, 0, 0, 0; ...
                                    -1, 1, 0, 0, 0, 0, 0, 0, 0, 0; ...
                                     0, 0, 1,-1, 0, 0,-1, 1, 0, 0; ...
                                     0, 0, 0, 0,-1, 0, 1,-1,-1, 1; ...
                                     0, 0, 0, 0  0,-1, 0, 0, 1,-1; ...
                                     ];
ThreeCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on';
                                            'kd*RNA'; 'kd*D_RNA'; 'kd*D1_RNA'; 
                                            'D*RNA'; 'D*D_RNA'; 'D*D_RNA'; 'D*D1_RNA'};
ThreeCompartmentModel.initialCondition = [0;1;0;0;0];
ThreeCompartmentModel.summarizeModel
ThreeCompartmentModel = ThreeCompartmentModel.formPropensitiesGeneral('TwoComparment_SpatialModel');
ThreeCompartmentModel.tSpan = linspace(0,10,21);
ThreeCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on','gene_off'};
save('ThreeComparment_SpatialModel','ThreeCompartmentModel')


% Four Compartment Model
FourCompartmentModel = SSIT();
FourCompartmentModel.parameters = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
FourCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'; 'D1_RNA'; 'D2_RNA'};
FourCompartmentModel.stoichiometry = [1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
                                     -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
                                      0, 0, 1,-1, 0, 0, 0,-1, 1, 0, 0, 0, 0; ...
                                      0, 0, 0, 0,-1, 0, 0, 1,-1,-1, 1, 0, 0; ...
                                      0, 0, 0, 0  0,-1, 0, 0, 0, 1,-1,-1, 1; ...
                                      0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1,-1
                                     ];
FourCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; ...
                                        'kd*RNA'; 'kd*D_RNA'; 'kd*D1_RNA'; 'kd*D2_RNA'; ...
                                        'D*RNA'; 'D*D_RNA'; 'D*D_RNA'; 'D*D1_RNA'; 'D*D1_RNA'; 'D*D2_RNA'};
FourCompartmentModel.initialCondition = [0;1;0;0;0;0;];
FourCompartmentModel.summarizeModel
FourCompartmentModel = FourCompartmentModel.formPropensitiesGeneral('TwoComparment_SpatialModel');
FourCompartmentModel.tSpan = linspace(0,10,21);
FourCompartmentModel.pdoOptions.unobservedSpecies = {'gene_on','gene_off'};
save('FourComparment_SpatialModel','FourCompartmentModel')


%% Calc Shit
num_compartments = [1, 2, 3, 4];
y1 = zeros([1, length(num_compartments)]);
y2 = zeros([1, length(num_compartments)]);

pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.makePlots = false;
pipelineArgs.pars = {'kon',1; 'koff',1; 'kr',100; 'kd',1};
saveName = 'OneCompartmentModel_regime1';
SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);
load([saveName, '.mat'])
y1(1) = outputs.determinates{1, "known_determinate"};
y2(1) = outputs.determinates{1, "unknown_determinate"};

pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = 4;
pipelineArgs.makePlots = false;
pipelineArgs.pars = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
saveName = 'TwoCompartmentModel_regime1';
SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);
load([saveName, '.mat'])
y1(2) = outputs.determinates{1, "known_determinate"};
y2(2) = outputs.determinates{1, "unknown_determinate"};

pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = 4;
pipelineArgs.makePlots = false;
pipelineArgs.pars = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
saveName = 'ThreeCompartmentModel_regime1';
SSIT('ThreeComparment_SpatialModel','ThreeCompartmentModel',{},pipeline,pipelineArgs,saveName);
load([saveName, '.mat'])
y1(3) = outputs.determinates{1, "known_determinate"};
y2(3) = outputs.determinates{1, "unknown_determinate"};

pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = 4;
pipelineArgs.makePlots = false;
pipelineArgs.pars = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
saveName = 'FourCompartmentModel_regime1';
SSIT('FourComparment_SpatialModel','FourCompartmentModel',{},pipeline,pipelineArgs,saveName);
load([saveName, '.mat'])
y1(4) = outputs.determinates{1, "known_determinate"};
y2(4) = outputs.determinates{1, "unknown_determinate"};






















