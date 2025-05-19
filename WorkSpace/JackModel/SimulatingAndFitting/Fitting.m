%% Fitting Simulation to multistate model
clear
close all 
addpath(genpath('../../src'));

%% Build Models 
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


% TODO change this to a different sized model
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




%% Load and parse data








%% Fit





















