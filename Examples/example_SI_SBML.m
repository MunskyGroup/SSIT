%% SSIT/Examples/example_SI_SBML

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SBML: Systems Biology Markup Language
%   * Export an SSIT model to an SBML file
%   * Create an SSIT model by loading an SBML file
%   (Requires the MATLAB SimBiology toolbox)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the Bursting Gene model from example_01_CreateSSITModels
% clear
% close all
% addpath(genpath('../src'));

% example_01_CreateSSITModels

% Load the Bursting Gene model created in example_01_CreateSSITModels:
load('ExampleSaveFiles/example_1_CreateSSITModels.mat','Model')

% View model summary:
Model.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Export the Bursting Gene model to an SBML file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of the SBML file to create:
sbmlFile = 'ExampleSaveFiles/BurstingGene_SBML.xml';

% Convert the model to SimBiology and write it to the SBML file:
Model.exportToSBML(sbmlFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Create a new SSIT model by loading the SBML file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create an empty SSIT model and load the SBML file into it:
Model_SBML = SSIT;
Model_SBML = Model_SBML.createModelFromSBML(sbmlFile);

% View model summary (should match the original Bursting Gene model):
Model_SBML.summarizeModel

% Set the times at which distributions will be computed:
Model_SBML.tSpan = linspace(0,50,101);

% Compile and store the reaction propensities:
Model_SBML = Model_SBML.formPropensitiesGeneral('Model_SBML');

% Solve with FSP:
Model_SBML = Model_SBML.solve(solver='FSP');

% Plot means and standard deviations:
Model_SBML.plotFSP(speciesNames=Model_SBML.species,...
    plotType='meansAndDevs', lineProps={'linewidth',4},...
    Title='Bursting Gene (from SBML)', TitleFontSize=26,...
    AxisLabelSize=20, TickLabelSize=20, YLabel='Molecule Count',...
    XLabel='Time', LegendFontSize=20, LegendLocation='southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Load a model from the SBML test suite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The second argument (scaleVolume=true) rescales non-integer initial
% amounts to integer molecule counts:
Model_Test = SSIT();
Model_Test = Model_Test.createModelFromSBML(...
    '../SBML_test_cases/00010/00010-sbml-l1v2.xml',true);
Model_Test = Model_Test.formPropensitiesGeneral('SBMLModel');
Model_Test.summarizeModel;
Model_Test = Model_Test.solve(solver='FSP');
Model_Test.plotFSP(speciesNames=Model_Test.species, plotType='meansAndDevs')
