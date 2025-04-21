%% Generates D vs |FIM| plot
clear
close all 
addpath(genpath('../../src'));
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
save('OneComparment_SpatialModel','OneCompartmentModel')

% Two Compartment Model
TwoCompartmentModel = SSIT();
TwoCompartmentModel.parameters = {'kon',1;'koff',1;'kr',100; 'D', 1; 'kd', 1};
TwoCompartmentModel.species = {'gene_on';'gene_off'; 'RNA'; 'D_RNA'};
TwoCompartmentModel.stoichiometry = [1, -1, 0, 0, 0, 0;... % species x reaction
                                    -1, 1, 0, 0, 0, 0;...
                                    0, 0, 1, -1, 0, -1;...
                                    0, 0, 0, 1, -1, 0;];
TwoCompartmentModel.propensityFunctions = {'kon*gene_off';'koff*gene_on';'kr*gene_on'; 'D*RNA'; 'kd*D_RNA'; 'kd*RNA'};
TwoCompartmentModel.initialCondition = [0;1;0;0;];
TwoCompartmentModel.summarizeModel
TwoCompartmentModel = TwoCompartmentModel.formPropensitiesGeneral('TwoComparment_SpatialModel');
TwoCompartmentModel.tSpan = linspace(0,100,21);
save('TwoComparment_SpatialModel','TwoCompartmentModel')


%% Pipeline
% specify fitting routine
pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = TwoCompartmentModel.parameters;
pipelineArgs.makePlots = true;

% Specify save name
saveName = 'GRModelFitResults';
model = SSIT('TwoComparment_SpatialModel','TwoCompartmentModel',{},pipeline,pipelineArgs,saveName);


% specify fitting routine
pipeline = 'multiModelFIMPipeline';
pipelineArgs.param_of_interest_index = NaN;
pipelineArgs.pars = OneCompartmentModel.parameters;
pipelineArgs.makePlots = true;

% Specify save name
saveName = 'GRModelFitResults';
model = SSIT('OneComparment_SpatialModel','OneCompartmentModel',{},pipeline,pipelineArgs,saveName);



%% Run Multiple Parameters 





























